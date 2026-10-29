/* Stream native SRA spots; reuse the accepted chrM k-mer candidate definition. */
#include <vdb/manager.h>
#include <vdb/database.h>
#include <vdb/table.h>
#include <vdb/cursor.h>
#include <zlib.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fs=std::filesystem;
void ok(rc_t rc,const char *where) { if(rc) throw std::runtime_error(std::string(where)+" rc="+std::to_string(rc)); }
template<typename T> size_t value(const VCursor* c,int64_t row,uint32_t col,T* dst,size_t cap) {
 uint32_t bits=0,offset=0,n=0;const void *p=nullptr;
 ok(VCursorCellDataDirect(c,row,col,&bits,&p,&offset,&n),"VCursorCellDataDirect");
 if(bits!=8*sizeof(T)||offset!=0||n>cap)throw std::runtime_error("unexpected VDB value representation/size");
 if(n)std::memcpy(dst,p,n*sizeof(T));
 return n;
}
std::vector<uint64_t> load_kmers(const std::string &path) {
 std::ifstream in(path,std::ios::binary);uint64_t n=0;in.read(reinterpret_cast<char*>(&n),8);
 if(!in||n==0||n>100000)throw std::runtime_error("invalid kmer header");
 std::vector<uint64_t> v(n);in.read(reinterpret_cast<char*>(v.data()),n*8);
 if(!in||in.peek()!=EOF||!std::is_sorted(v.begin(),v.end()))throw std::runtime_error("invalid kmer set");
 return v;
}
bool encode21(const char *s,uint64_t &v) {
 v=0;
 for(int i=0;i<21;i++) {
  unsigned char ch=std::toupper(static_cast<unsigned char>(s[i]));
  int x=ch=='A'?0:ch=='C'?1:ch=='G'?2:ch=='T'?3:ch=='N'?4:-1;
  if(x<0)return false;
  v=(v<<3)|static_cast<uint64_t>(x);
 }
 return true;
}
bool mitochondrial_candidate(const char *seq,size_t n,const std::vector<uint64_t> &kmers) {
 if(n<50)return false;
 for(size_t j=0;j+21<=n;j+=10) {
  uint64_t v;if(encode21(seq+j,v)&&std::binary_search(kmers.begin(),kmers.end(),v))return true;
 }
 uint64_t v;return encode21(seq+n-21,v)&&std::binary_search(kmers.begin(),kmers.end(),v);
}
std::unordered_set<std::string> load_axis(const std::string &path) {
 gzFile f=gzopen(path.c_str(),"rb");if(!f)throw std::runtime_error("cannot open barcode gzip");
 std::unordered_set<std::string> axis;char line[128];
 while(gzgets(f,line,sizeof line)) {
  std::string s(line);while(!s.empty()&&(s.back()=='\n'||s.back()=='\r'))s.pop_back();
  if(s.size()!=18||s.substr(16)!="-1")throw std::runtime_error("invalid called barcode");
  if(!axis.insert(s.substr(0,16)).second)throw std::runtime_error("duplicate called barcode");
 }
 if(gzclose(f)!=Z_OK||axis.size()!=4379)throw std::runtime_error("called barcode gzip mismatch");
 return axis;
}
std::unordered_map<std::string,char> load_roles(const std::string &path,const std::unordered_set<std::string>&axis) {
 std::ifstream in(path);std::string line;std::getline(in,line);
 if(line!="barcode\tcluster\trole")throw std::runtime_error("invalid role header");
 std::unordered_map<std::string,char> roles;int sources=0,targets=0;
 while(std::getline(in,line)) {
  std::istringstream row(line);std::string cell,cluster,role;
  std::getline(row,cell,'\t');std::getline(row,cluster,'\t');std::getline(row,role);
  if(cell.size()!=18||cell.substr(16)!="-1")throw std::runtime_error("invalid role barcode");
  cell=cell.substr(0,16);char kind=role=="source"?'S':role=="target"?'T':'?';
  if(!axis.count(cell)||kind=='?'||!roles.emplace(cell,kind).second)throw std::runtime_error("invalid/duplicate role");
  sources+=kind=='S';targets+=kind=='T';
 }
 if(sources!=108||targets!=2001)throw std::runtime_error("fixed role count mismatch");
 return roles;
}
struct Counts {
 uint64_t spots=0,valid_structure=0,unexpected_structure=0,barcode_acgt=0,axis_match=0,
          barcode_q20=0,role_match=0,source=0,target=0,umi_acgt=0,umi_q20=0,candidates=0;
};
void add(Counts &a,const Counts &b) {
 a.spots+=b.spots;a.valid_structure+=b.valid_structure;a.unexpected_structure+=b.unexpected_structure;
 a.barcode_acgt+=b.barcode_acgt;a.axis_match+=b.axis_match;a.barcode_q20+=b.barcode_q20;
 a.role_match+=b.role_match;a.source+=b.source;a.target+=b.target;a.umi_acgt+=b.umi_acgt;
 a.umi_q20+=b.umi_q20;a.candidates+=b.candidates;
}
bool valid_acgt(const char *p,size_t n) {
 for(size_t i=0;i<n;i++)if(p[i]!='A'&&p[i]!='C'&&p[i]!='G'&&p[i]!='T')return false;
 return true;
}
int main(int argc,char **argv) try {
 if(argc==3&&std::string(argv[1])=="--kmer-test") {
  auto kmers=load_kmers(argv[2]);std::string s;while(std::getline(std::cin,s))std::cout<<mitochondrial_candidate(s.data(),s.size(),kmers)<<"\n";return 0;
 }
 if(argc!=10)throw std::runtime_error("args: archive run expected_spots offset axis.gz roles.tsv kmers.bin output_dir batch_spots");
 std::string archive=argv[1],run=argv[2],outdir=argv[8];uint64_t expected=std::stoull(argv[3]),offset=std::stoull(argv[4]),batch_size=std::stoull(argv[9]);
 if(!batch_size)throw std::runtime_error("batch size zero");
 auto kmers=load_kmers(argv[7]);auto axis=load_axis(argv[5]);auto roles=load_roles(argv[6],axis);
 fs::create_directories(outdir);fs::path ledger=fs::path(outdir)/(run+"_batches.csv");
 if(fs::exists(ledger))throw std::runtime_error("existing batch ledger; preserve before any resume");
 std::ofstream log(ledger);if(!log)throw std::runtime_error("cannot create ledger");
 log<<"run,first_spot,last_spot,spots,valid_structure,unexpected_structure,barcode_acgt,axis_match,barcode_q20,role_match,source,target,umi_acgt,umi_q20,candidates,candidate_file,seconds\n";
 const VDBManager *mgr=nullptr;const VDatabase *db=nullptr;const VTable *table=nullptr;const VCursor *cursor=nullptr;
 ok(VDBManagerMakeRead(&mgr,nullptr),"manager");ok(VDBManagerOpenDBRead(mgr,&db,nullptr,"%s",archive.c_str()),"open archive");
 ok(VDatabaseOpenTableRead(db,&table,"SEQUENCE"),"open sequence");ok(VTableCreateCursorRead(table,&cursor),"cursor");
 uint32_t ir,iq,il,it,is;
 ok(VCursorAddColumn(cursor,&ir,"READ"),"READ");ok(VCursorAddColumn(cursor,&iq,"QUALITY"),"QUALITY");
 ok(VCursorAddColumn(cursor,&il,"READ_LEN"),"READ_LEN");ok(VCursorAddColumn(cursor,&it,"READ_TYPE"),"READ_TYPE");ok(VCursorAddColumn(cursor,&is,"READ_START"),"READ_START");
 ok(VCursorOpen(cursor),"open cursor");int64_t first;uint64_t count;ok(VCursorIdRange(cursor,ir,&first,&count),"id range");
 if(first!=1||count!=expected)throw std::runtime_error("whole-run spot count mismatch");
 Counts batch,total;gzFile gz=nullptr;std::string temp,final;uint64_t batch_first=1,batch_no=0;auto started=std::chrono::steady_clock::now();auto batch_time=started;
 for(uint64_t i=0;i<count;i++) {
  int64_t row=first+static_cast<int64_t>(i);char seq[512];uint8_t quality[512],types[16];uint32_t lens[16],starts[16];
  size_t nr=value(cursor,row,ir,seq,512),nq=value(cursor,row,iq,quality,512),nl=value(cursor,row,il,lens,16),nt=value(cursor,row,it,types,16),ns=value(cursor,row,is,starts,16);
  batch.spots++;
  bool good=nl==4&&nt==4&&ns==4&&nr==139&&nq==139&&
            lens[0]==10&&lens[1]==10&&lens[2]==28&&lens[3]==91&&starts[0]==0&&starts[1]==10&&starts[2]==20&&starts[3]==48&&
            !(types[0]&1)&&!(types[1]&1)&&!(types[2]&1)&&(types[3]&1);
  if(!good)batch.unexpected_structure++;
  else {
   batch.valid_structure++;char *barcode=seq+20,*umi=seq+36,*rna=seq+48;
   if(valid_acgt(barcode,16)) {
    batch.barcode_acgt++;std::string cell(barcode,16);
    if(axis.count(cell)) {
     batch.axis_match++;bool bq=true;for(int z=20;z<36;z++)bq&=quality[z]>=20;
     if(bq)batch.barcode_q20++;
     auto it=roles.find(cell);
     if(it!=roles.end()) {
      batch.role_match++;if(it->second=='S')batch.source++;else batch.target++;
      if(valid_acgt(umi,12)) {
       batch.umi_acgt++;bool uq=true;for(int z=36;z<48;z++)uq&=quality[z]>=20;
       if(uq) {
        batch.umi_q20++;
        if(mitochondrial_candidate(rna,91,kmers)) {
         batch.candidates++;
         if(!gz) {
          std::ostringstream name;name<<run<<"_batch_"<<batch_no<<".fastq.gz";
          final=(fs::path(outdir)/name.str()).string();temp=final+".writing";
          if(fs::exists(final)||fs::exists(temp))throw std::runtime_error("candidate batch exists");
          gz=gzopen(temp.c_str(),"wb2");if(!gz)throw std::runtime_error("cannot create candidate gzip");
         }
         char q[92];for(int z=0;z<91;z++) { if(quality[48+z]>93)throw std::runtime_error("invalid phred");q[z]=static_cast<char>(quality[48+z]+33); }q[91]=0;
         std::string cdna(rna,91);std::transform(cdna.begin(),cdna.end(),cdna.begin(),[](unsigned char x){return std::toupper(x);});
         auto id=offset+static_cast<uint64_t>(row);
         if(gzprintf(gz,"@%s_%.*s_%llu\n%s\n+\n%s\n",cell.c_str(),12,umi,(unsigned long long)id,cdna.c_str(),q)<=0)
          throw std::runtime_error("candidate gzip write failed");
        }
       }
      }
     }
    }
   }
  }
  if((i+1)%1000000==0)std::cerr<<run<<" spots="<<(i+1)<<" candidates="<<(total.candidates+batch.candidates)<<"\n";
  if((i+1)%batch_size==0||i+1==count) {
   if(gz) { if(gzclose(gz)!=Z_OK)throw std::runtime_error("gzip close failed");gz=nullptr;fs::rename(temp,final); }
   double secs=std::chrono::duration<double>(std::chrono::steady_clock::now()-batch_time).count();
   log<<run<<","<<batch_first<<","<<(i+1)<<","<<batch.spots<<","<<batch.valid_structure<<","<<batch.unexpected_structure<<","<<batch.barcode_acgt<<","<<batch.axis_match<<","<<batch.barcode_q20<<","<<batch.role_match<<","<<batch.source<<","<<batch.target<<","<<batch.umi_acgt<<","<<batch.umi_q20<<","<<batch.candidates<<","<<(batch.candidates?fs::path(final).filename().string():"")<<","<<secs<<"\n";
   log.flush();if(!log)throw std::runtime_error("batch ledger write failed");add(total,batch);batch=Counts{};batch_first=i+2;batch_no++;batch_time=std::chrono::steady_clock::now();final.clear();temp.clear();
  }
 }
 if(total.spots!=count||total.valid_structure+total.unexpected_structure!=count)throw std::runtime_error("spot accounting mismatch");
 std::cerr<<run<<" COMPLETE spots="<<total.spots<<" candidates="<<total.candidates<<" unexpected="<<total.unexpected_structure
          <<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"\n";
 VCursorRelease(cursor);VTableRelease(table);VDatabaseRelease(db);VDBManagerRelease(mgr);
 return 0;
} catch(const std::exception &e){std::cerr<<"ERROR "<<e.what()<<"\n";return 1;}
