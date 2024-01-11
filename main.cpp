//----------------------------------------------------------------
// Name        : main.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description :
//----------------------------------------------------------------
#include <iostream>
#include <map>
#include <getopt.h>
#include "htslib/sam.h"
//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------
#define bit_set(w,m)    ((w) |= (m))
#define bit_res(w,m)    ((w) &= ~(m))
#define bit_test(w,m)   ((w) & (m))
//----------------------------------------------------------------
#define INPUT_BAM_FILE      'i'
#define OUTPUT_BAM_FILE     'o'
#define MAX_ALIGNMENT_SCORE 's'
#define HARD_CLIPPING       'c'
#define HELP                'h'
#define SHORT_OPTIONS       "i:o:s:ch"
//----------------------------------------------------------------
struct option longOptions[] =
{
    {"input-bam-file",required_argument,nullptr,INPUT_BAM_FILE},
    {"output-bam-file",required_argument,nullptr,OUTPUT_BAM_FILE},
    {"max-alignment-score",required_argument,nullptr,MAX_ALIGNMENT_SCORE},
    {"hard-clipping",no_argument,nullptr,HARD_CLIPPING},
    {"help",no_argument,nullptr,HELP}
};
//----------------------------------------------------------------
bool SplitSingleRead(const bam1_t * read,multimap<hts_pos_t,bam1_t*> & readBuffer,bool useHardClipping)
{
    bam1_t * newRead;

    const char * queryName=bam_get_qname(read);
    size_t queryNameLen=strlen(bam_get_qname(read));

    uint16_t flag;

    hts_pos_t startPos,pos; startPos=pos=read->core.pos;

    size_t newCigarLen=0;
    uint32_t newCigar[256];

    size_t newQueryLen=0;
    char newSeq[16384];
    char newQual[16384];

    uint32_t cigarLen=read->core.n_cigar;
    uint32_t * cigar=bam_get_cigar(read);

    size_t queryIndex=0;
    uint8_t * seq=(uint8_t*)bam_get_seq(read);
    char * qual=(char*)bam_get_qual(read);

    size_t l_aux=bam_get_l_aux(read);
    uint8_t * aux=bam_get_aux(read);

    for(uint32_t i=0;i<cigarLen;i++)
    {
        uint32_t op=bam_cigar_op(cigar[i]);
        uint32_t opLen=uint32_t(bam_cigar_oplen(cigar[i]));

        switch(op)
        {
            //----------------------------------------------------------------
            //Alignment match and mismatch consumes reference and query
            //----------------------------------------------------------------

            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:

                newCigar[newCigarLen++]=bam_cigar_gen(opLen,op);

                pos+=opLen;

                for(;opLen>0;queryIndex++,newQueryLen++,opLen--)
                {
                    newSeq[newQueryLen] = seq_nt16_str[bam_seqi(seq,queryIndex)];
                    newQual[newQueryLen] = qual[queryIndex];
                }

                continue;

            //----------------------------------------------------------------
            //Insertion and soft clipping only consume query
            //----------------------------------------------------------------

            case BAM_CINS:
            case BAM_CSOFT_CLIP:

                newCigar[newCigarLen++]=bam_cigar_gen(opLen,op);

                for(;opLen>0;queryIndex++,newQueryLen++,opLen--)
                {
                    newSeq[newQueryLen] = seq_nt16_str[bam_seqi(seq,queryIndex)];
                    newQual[newQueryLen] = qual[queryIndex];
                }

                continue;

            //----------------------------------------------------------------
            //Deletion consumes reference but not query
            //----------------------------------------------------------------

            case BAM_CDEL:

                newCigar[newCigarLen++]=bam_cigar_gen(opLen,op);
                pos+=opLen;
                continue;

            //----------------------------------------------------------------
            //Hard clipping and padding consumes nothing simply copy the cigar op
            //----------------------------------------------------------------

            case BAM_CHARD_CLIP:
            case BAM_CPAD:

                newCigar[newCigarLen++]=bam_cigar_gen(opLen,op);
                continue;

            //----------------------------------------------------------------
            //Split on skipped region from the reference (N)
            //----------------------------------------------------------------

            case BAM_CREF_SKIP:

                if(useHardClipping==true) newCigar[newCigarLen++]=bam_cigar_gen(opLen-(opLen>>1),BAM_CHARD_CLIP);

                flag=read->core.flag;
                bit_set(flag,BAM_FPAIRED);
                bit_set(flag,BAM_FPROPER_PAIR);
                bit_res(flag,BAM_FMUNMAP);
                bit_test(flag,BAM_FREVERSE)!=0 ? bit_set(flag,BAM_FMREVERSE) : bit_res(flag,BAM_FMREVERSE);
                bit_set(flag,BAM_FREAD1);
                if(startPos!=read->core.pos){bit_set(flag,BAM_FREAD2);} else{bit_res(flag,BAM_FREAD2);}

                pos+=opLen;

                newRead=bam_init1();

                if(bam_set1(newRead,queryNameLen,queryName,flag,read->core.tid,startPos,read->core.qual,newCigarLen,newCigar,read->core.tid,pos,read->core.isize,newQueryLen,newSeq,newQual,l_aux)<0)
                {
                    cerr << "Error: Could not initialize new read" << endl;
                    bam_destroy1(newRead);
                    return false;
                }

                newRead->l_data+=l_aux;
                memcpy(bam_get_aux(newRead),aux,l_aux);

                readBuffer.insert({newRead->core.pos,newRead});  //Store the read in the read buffer

                startPos=pos;

                newCigarLen=0;
                newQueryLen=0;

                if(useHardClipping==true && opLen>1) newCigar[newCigarLen++]=bam_cigar_gen(opLen>>1,BAM_CHARD_CLIP);

                continue;
        }
    }

    flag=read->core.flag;

    if(startPos!=read->core.pos)
    {
        bit_set(flag,BAM_FPAIRED);
        bit_set(flag,BAM_FPROPER_PAIR);
        bit_res(flag,BAM_FMUNMAP);
        if(bit_test(flag,BAM_FREVERSE)!=0){bit_set(flag,BAM_FMREVERSE);} else{bit_res(flag,BAM_FMREVERSE);}
        bit_res(flag,BAM_FREAD1);
        bit_set(flag,BAM_FREAD2);
    }

    newRead=bam_init1();

    if(bam_set1(newRead,queryNameLen,queryName,flag,read->core.tid,startPos,read->core.qual,newCigarLen,newCigar,read->core.mtid,read->core.mpos,read->core.isize,newQueryLen,newSeq,newQual,l_aux)<0)
    {
        cerr << "Error: Could not initialize new read" << endl;
        bam_destroy1(newRead);
        return false;
    }

    newRead->l_data+=l_aux;
    memcpy(bam_get_aux(newRead),aux,l_aux);

    readBuffer.insert({newRead->core.pos,newRead});  //Store the read in the read buffer (re-order buffer)

    return true;
}
//----------------------------------------------------------------
bool SplitMultiRead(const bam1_t * read,multimap<hts_pos_t,bam1_t*> & readBuffer,bool useHardClipping)
{
    bam1_t * newRead;

    const char * queryName=bam_get_qname(read);
    size_t queryNameLen=strlen(bam_get_qname(read));

    uint16_t flag,newFlag;

    hts_pos_t startPos,pos; startPos=pos=read->core.pos;

    size_t newCigarLen=0;
    uint32_t newCigar[256];

    size_t newQueryLen=0;
    char newSeq[16384];
    char newQual[16384];

    uint32_t cigarLen=read->core.n_cigar;
    uint32_t * cigar=bam_get_cigar(read);

    size_t queryIndex=0;
    uint8_t * seq=(uint8_t*)bam_get_seq(read);
    char * qual=(char*)bam_get_qual(read);

    size_t l_aux=bam_get_l_aux(read);
    uint8_t * aux=bam_get_aux(read);

    for(uint32_t i=0;i<cigarLen;i++)
    {
        uint32_t op=bam_cigar_op(cigar[i]);
        uint32_t opLen=uint32_t(bam_cigar_oplen(cigar[i]));

        switch(op)
        {
            //----------------------------------------------------------------
            //Alignment match and mismatch consumes reference and query
            //----------------------------------------------------------------

            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:

                newCigar[newCigarLen++]=bam_cigar_gen(opLen,op);

                pos+=opLen;

                for(;opLen>0;queryIndex++,newQueryLen++,opLen--)
                {
                    newSeq[newQueryLen] = seq_nt16_str[bam_seqi(seq,queryIndex)];
                    newQual[newQueryLen] = qual[queryIndex];
                }

                continue;

            //----------------------------------------------------------------
            //Insertion and soft clipping only consume query
            //----------------------------------------------------------------

            case BAM_CINS:
            case BAM_CSOFT_CLIP:

                newCigar[newCigarLen++]=bam_cigar_gen(opLen,op);

                for(;opLen>0;queryIndex++,newQueryLen++,opLen--)
                {
                    newSeq[newQueryLen] = seq_nt16_str[bam_seqi(seq,queryIndex)];
                    newQual[newQueryLen] = qual[queryIndex];
                }

                continue;

            //----------------------------------------------------------------
            //Deletion consumes reference but not query
            //----------------------------------------------------------------

            case BAM_CDEL:

                newCigar[newCigarLen++]=bam_cigar_gen(opLen,op);
                pos+=opLen;
                continue;

            //----------------------------------------------------------------
            //Hard clipping and padding consumes nothing simply copy the cigar op
            //----------------------------------------------------------------

            case BAM_CHARD_CLIP:
            case BAM_CPAD:

                newCigar[newCigarLen++]=bam_cigar_gen(opLen,op);
                continue;

            //----------------------------------------------------------------
            //Split on skipped region from the reference (N)
            //----------------------------------------------------------------

            case BAM_CREF_SKIP:

                if(useHardClipping==true) newCigar[newCigarLen++]=bam_cigar_gen(opLen-(opLen>>1),BAM_CHARD_CLIP);

                flag=newFlag=read->core.flag;
                bit_res(newFlag,BAM_FMUNMAP);
                if(bit_test(flag,BAM_FREVERSE)!=0) {bit_set(newFlag,BAM_FMREVERSE);} else{ bit_res(newFlag,BAM_FMREVERSE);}
                if(bit_test(flag,BAM_FREAD2)!=0) bit_set(newFlag,BAM_FREAD1);
                if(startPos!=read->core.pos && bit_test(flag,BAM_FREAD1)!=0) bit_set(newFlag,BAM_FREAD2);

                pos+=opLen;

                newRead=bam_init1();

                if(bam_set1(newRead,queryNameLen,queryName,newFlag,read->core.tid,startPos,read->core.qual,newCigarLen,newCigar,read->core.tid,pos,read->core.isize,newQueryLen,newSeq,newQual,l_aux)<0)
                {
                    cerr << "Error: Could not initialize new read" << endl;
                    bam_destroy1(newRead);
                    return false;
                }

                newRead->l_data+=l_aux;
                memcpy(bam_get_aux(newRead),aux,l_aux);

                readBuffer.insert({newRead->core.pos,newRead});

                startPos=pos;

                newCigarLen=0;
                newQueryLen=0;

                if(useHardClipping==true  && opLen>1) newCigar[newCigarLen++]=bam_cigar_gen(opLen>>1,BAM_CHARD_CLIP);

                continue;
        }
    }

    newFlag=read->core.flag;

    if(startPos!=read->core.pos)
    {
        if(bit_test(newFlag,BAM_FREAD1)!=0) bit_set(newFlag,BAM_FREAD2);
    }

    newRead=bam_init1();

    if(bam_set1(newRead,queryNameLen,queryName,newFlag,read->core.tid,startPos,read->core.qual,newCigarLen,newCigar,read->core.mtid,read->core.mpos,read->core.isize,newQueryLen,newSeq,newQual,l_aux)<0)
    {
        cerr << "Error: Could not initialize new read" << endl;
        bam_destroy1(newRead);
        return false;
    }

    newRead->l_data+=l_aux;
    memcpy(bam_get_aux(newRead),aux,l_aux);

    readBuffer.insert({newRead->core.pos,newRead});

    return true;
}
//----------------------------------------------------------------
bool FlushReadBuffer(samFile * bamFile,bam_hdr_t * bamHeader,multimap<hts_pos_t,bam1_t*> & readBuffer)
{
    auto itEnd=readBuffer.end();

    for(auto it=readBuffer.begin();it!=itEnd;it++)
    {
        if(sam_write1(bamFile,bamHeader,it->second)<0)
        {
            cerr << "Error: Could not write alignment" << endl;
            for(;it!=itEnd;it++){bam_destroy1(it->second);} readBuffer.clear();
            return false;
        }

        bam_destroy1(it->second);
    }

    readBuffer.clear();

    return true;
}
//----------------------------------------------------------------
int Run(const string & inputBamFileName,const string & outputBamFileName,uint8_t maxAlignmentScore,bool useHardClipping)
{
    int ret=1;

    samFile * inputBamFile=nullptr;
    bam_hdr_t * bamHeader=nullptr;
    samFile * outputBamFile=nullptr;

    int32_t prevTID=-1;
    hts_pos_t prevPos=-1;

    bam1_t * read=nullptr;
    multimap<hts_pos_t,bam1_t*> readBuffer;

    //----------------------------------------------------------------
    //Open input bam file
    //----------------------------------------------------------------

    inputBamFile=hts_open(inputBamFileName.c_str(),"r");

    if(inputBamFile==nullptr)
    {
        cerr << "Error: Could not open input bam file" << endl;
        goto CLEAN_UP;
    }

    //----------------------------------------------------------------
    //Read bam file header
    //----------------------------------------------------------------

    bamHeader=sam_hdr_read(inputBamFile);

    if(bamHeader==nullptr)
    {
        cerr << "Error: Could not read bam header" << endl;
        goto CLEAN_UP;
    }

    //----------------------------------------------------------------
    //Create output bam file
    //----------------------------------------------------------------

    outputBamFile=hts_open(outputBamFileName.c_str(),"wb");

    if(outputBamFile==nullptr)
    {
        cerr << "Error: Could not create output bam file" << endl;
        goto CLEAN_UP;
    }

    //----------------------------------------------------------------
    //Write bam file header
    //----------------------------------------------------------------

    if(sam_hdr_write(outputBamFile,bamHeader)<0)
    {
        cerr << "Error: Could not write bam header" << endl;
        goto CLEAN_UP;
    }

    //----------------------------------------------------------------
    //Process files
    //----------------------------------------------------------------

    read=bam_init1();

    while(sam_read1(inputBamFile,bamHeader,read)>=0)
    {
        //----------------------------------------------------------------
        //Adjust max alignment score (Do that for any read)
        //----------------------------------------------------------------

        if(read->core.qual>maxAlignmentScore) read->core.qual=maxAlignmentScore;

        //----------------------------------------------------------------
        //Simply copy unmapped reads (Flush the read buffer prior to writing out unmapped reads, they should all come together but the code does not enforce this)
        //----------------------------------------------------------------

        if(bit_test(read->core.flag,BAM_FUNMAP)!=0) //Unmapped reads are a special case
        {
            prevTID=-1;
            prevPos=-1;

            if(readBuffer.empty()==false)
            {
                if(FlushReadBuffer(outputBamFile,bamHeader,readBuffer)==false) goto CLEAN_UP;
            }

            if(sam_write1(outputBamFile,bamHeader,read)<0)
            {
                cerr << "Error: Could not write unmapped alignment" << endl;
                goto CLEAN_UP;
            }

            continue;
        }

        //----------------------------------------------------------------
        //Check if bam file is sorted on coordinate
        //----------------------------------------------------------------

        if((read->core.tid<prevTID) || (read->core.tid==prevTID && read->core.pos<prevPos))
        {
            cerr << "Error: Bam file is not sorted on coordinate (Please sort the bam file first)" << endl;
            goto CLEAN_UP;
        }

        //----------------------------------------------------------------
        //Change of chromosome? Flush the buffer with reads of previous chromosome
        //----------------------------------------------------------------

        if(read->core.tid>prevTID)
        {
            prevTID=read->core.tid;

            if(readBuffer.empty()==false)
            {
                if(FlushReadBuffer(outputBamFile,bamHeader,readBuffer)==false) goto CLEAN_UP;   //Code path read buffer clear
            }
        }

        prevPos=read->core.pos;

        //----------------------------------------------------------------
        //Template with single segment
        //----------------------------------------------------------------

        if(bit_test(read->core.flag,BAM_FPAIRED)==0)
        {
            if(SplitSingleRead(read,readBuffer,useHardClipping)==false) goto CLEAN_UP;
        }

        //----------------------------------------------------------------
        //Template with multiple segment
        //----------------------------------------------------------------

        else
        {
            if(SplitMultiRead(read,readBuffer,useHardClipping)==false) goto CLEAN_UP;
        }

        //----------------------------------------------------------------
        //Syncrhonize the read buffer to disk up until prevPos (Special case of FlushReadBuffer)
        //----------------------------------------------------------------

        auto itBegin=readBuffer.begin();
        auto itEnd=readBuffer.upper_bound(prevPos);

        for(auto it=itBegin;it!=itEnd;it++)
        {
            if(sam_write1(outputBamFile,bamHeader,it->second)<0)
            {
                cerr << "Error: Could not write read" << endl;
                for(itEnd=readBuffer.end();it!=itEnd;it++){bam_destroy1(it->second);} readBuffer.clear();   //Writing failed clean up the rest of the buffer
                goto CLEAN_UP;
            }

            bam_destroy1(it->second);
        }

        readBuffer.erase(itBegin,itEnd);
    }

    //----------------------------------------------------------------
    //Write out the last reads in the buffer if any
    //----------------------------------------------------------------

    if(readBuffer.empty()==false)
    {
        if(FlushReadBuffer(outputBamFile,bamHeader,readBuffer)==false) goto CLEAN_UP;
    }

    ret=0;

    //----------------------------------------------------------------
    //Clean up
    //----------------------------------------------------------------

CLEAN_UP:;

    for(const auto & it : readBuffer)bam_destroy1(it.second);     //No checks anymore on readBuffer no need to call clear
    if(read!=nullptr) bam_destroy1(read);
    if(outputBamFile!=nullptr) hts_close(outputBamFile);
    if(bamHeader!=nullptr) sam_hdr_destroy(bamHeader);
    if(inputBamFile!=nullptr) hts_close(inputBamFile);

    if(sam_index_build(outputBamFileName.c_str(),0)!=0) //Build index only when file is closed and finished
    {
        cerr << "Error: Could not create index" << endl;
        return 1;
    }

    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------

    return ret;
}
//----------------------------------------------------------------
int main(int argc, char** argv)
{
    //----------------------------------------------------------------
    //Get input arguments
    //----------------------------------------------------------------

    bool showHelp=(argc==1);
    int option,optionIndex;

    string inputBamFileName,outputBamFileName;
    int32_t maxAlignmentScore=255;
    bool useHardClipping=false;

    while((option=getopt_long(argc,argv,SHORT_OPTIONS,longOptions,&optionIndex))>=0)
    {
        switch(option)
        {
            case INPUT_BAM_FILE:

                inputBamFileName=string(optarg);
                break;

            case OUTPUT_BAM_FILE:

                outputBamFileName=string(optarg);
                break;

            case MAX_ALIGNMENT_SCORE:

                maxAlignmentScore=atoi(optarg);
                break;

            case HARD_CLIPPING:

                useHardClipping=true;
                break;

            case HELP:

                showHelp=true;
                break;
        }
    }

    //----------------------------------------------------------------
    //Show help
    //----------------------------------------------------------------

    if(showHelp)
    {
        cerr << "split_spliced_reads [options]"                                                         << endl;
        cerr                                                                                            << endl;
        cerr << "-i --input-bam-file <text>     Single input bam file (required)"                       << endl;
        cerr << "-o --output-bam-file <text>    Single output bam file (optional default _split.bam)"   << endl;
        cerr << "-s --max-alignment-score <int> Maximum alignment score (optional default 255)"         << endl;
        cerr << "-c --hard-clipping <void>      When introduced use hard clipping"                      << endl;
        cerr << "-h --help <void>               This help"                                              << endl;
        cerr                                                                                            << endl;

        return 0;
    }

    //----------------------------------------------------------------
    //Check input arguments
    //----------------------------------------------------------------

    if(inputBamFileName.empty())
    {
        cerr << "Error: Please specify an input bam file!" << endl;
        return 1;
    }

    if(outputBamFileName.empty())
    {
        outputBamFileName=inputBamFileName.substr(0,inputBamFileName.find_last_of('.')) + string("_split.bam");
    }

    if(maxAlignmentScore<20) maxAlignmentScore=20;
    if(maxAlignmentScore>255) maxAlignmentScore=255;

    //----------------------------------------------------------------
    //Lets go
    //----------------------------------------------------------------

    return Run(inputBamFileName,outputBamFileName,uint8_t(maxAlignmentScore),useHardClipping);
}
//----------------------------------------------------------------


