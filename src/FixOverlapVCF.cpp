#include <iostream>
#include <vector>
#include <algorithm>
#include "htslib/hfile.h"
#include "htslib/htslib/vcf.h"
#include "htslib/htslib/hts.h"
#include <string>
#define NULL_int32_t -2147483648 // this is also min btw but hoe hts lib stores null
#define DEBUG false
#define DEBUG_V false


using std::cout;
using std::cerr;
using std::cin;

int print_help(void){
    cerr << "-v\t      \t If this tag is included send a message to /dev/stderr every 1k pos\n"
         << "-i\t<path>\t Path to input vcf file (it should be bgziped and indexed)\n"
         << "-o\t<path>\t Path to output fixed file to\n"
         << "-r\t<str> \t String for the reference tag in the vcf [RO]\n"
         << "-d\t<str> \t String for the depth tag in the vcf [DP]\n";
    exit(-1);
}


class args{
    public:
        char *                vcf_path  = nullptr;
        htsFile *             vcf_file  = nullptr;
        bcf_hdr_t *           hdr       = nullptr;
        htsFile *             out_vcf   = nullptr;
        std::string           ref_tag   = "RO";
        std::string           depth_tag = "DP";
        std::vector<int32_t>  depth ;
        std::vector<int32_t>  ref ;
        int32_t *             dp_dat    = nullptr;
        int32_t *             ref_dat   = nullptr;
        unsigned short int    verbose   = 0;
        int                   nsamps    = 0;
    int parse(int len , char ** input){
        int i = 1;
        int out_idx = 0; // we need this b/c we need to encuse we have the input for our header template
        while (i < len){
            // when we get the vcf file input we can initialise a bunch of stuff so we will do that here
            if (strncmp(input[i], "-i", 2) == 0){
                vcf_path = input[i + 1];
                vcf_file = bcf_open(input[i + 1], "r");
                if(vcf_file == nullptr){cerr << "Failed to open file " << vcf_path << std::endl; exit(-1);}
                hdr = bcf_hdr_read(vcf_file);
                if(hdr == nullptr){cerr << "Failed to open header of file " << vcf_path << std::endl; exit(-1);}
                nsamps = bcf_hdr_nsamples(hdr);
                depth.resize((unsigned int) nsamps);
                ref.resize((unsigned int) nsamps);
                i  += 2;
            } else if (strncmp(input[i], "-r", 2) == 0){
                ref_tag = input[i + 1];
                i += 2;
            } else if (strncmp(input[i], "-d", 2) == 0){
                depth_tag = input[i + 1];
                i += 2;
            } else if (strncmp(input[i], "-o", 2) == 0) {
                out_idx = i + 1;
                i += 2; // save for later
            } else if (strncmp(input[i], "-h", 2) == 0){
                print_help();
            } else if (strncmp(input[i], "-v", 2) == 0){
                verbose += 1000;
                i += 1;
            }else {
                cerr << "Unrecognised flag:\t" << input[i++] << std::endl;
            }

        }
        if(out_idx == 0 or hdr == nullptr){
            print_help();
        }
        out_vcf = bcf_open(input[out_idx],"w");
        if (bcf_hdr_write(out_vcf, hdr) != 0){cerr << "failed to copy header!!" << std::endl ; exit(-1);}
        dp_dat = depth.data();
        ref_dat = ref.data();
        return  0;
    }
};


void output_lines(args arg, std::vector<bcf1_t *>& ori_lines, std::vector<int>& max){
    int pos = max[0];
    int i = 0 ;
    while(i < ori_lines.size() ){
        // only print things if they don't overlap the next region...
        if(ori_lines[i]->pos + ori_lines[i]->rlen <= pos){
#if DEBUG_V
            cerr << "Did we update depth?: \t  " <<
#endif
            bcf_update_format_int32(arg.hdr, ori_lines[i], arg.depth_tag.c_str(), (void *) arg.depth.data(), arg.nsamps);
#if DEBUG_V
            cerr << "\nDid we update reference?: \t  " <<
#endif
            bcf_update_format_int32(arg.hdr, ori_lines[i], arg.ref_tag.c_str(), (void *) arg.ref.data(), arg.nsamps);
#if DEBUG_V
            cerr << "\n";
#endif
            if(bcf_write(arg.out_vcf, arg.hdr, ori_lines[i])) {
                cerr << "failed to write line to the bcf file\n";
                exit(-1);
            } else {
                // we have to get rid of anything we print here
                ori_lines.erase(ori_lines.begin() + i);
                max.erase(max.begin() + i);
            }
        } else {
            // this line overlaps, skip it !!
            i++;

        }
    }
    if(!ori_lines.empty()) {
#if DEBUG
        cerr << "using position: "  << ori_lines[0]->pos << "as a starting point\n";
#endif
        bcf_get_format_int32(arg.hdr, ori_lines[0], arg.depth_tag.c_str(), &arg.dp_dat, &arg.nsamps);
        bcf_get_format_int32(arg.hdr, ori_lines[0], arg.ref_tag.c_str(), &arg.ref_dat, &arg.nsamps);
        if (ori_lines.size() > 1) {
            cerr << "Theoretically you never get here...\n" ;
        }
    } else {
        // if we are not overwriting this we need to set it to null (Min value)
        std::fill(arg.depth.begin(), arg.depth.end(), INT32_MIN);
        std::fill(arg.ref.begin(), arg.ref.end(), INT32_MIN);
    }
}


int main(int argc, char ** argv) {
    args user_input;
    std::vector<bcf1_t *> lines ;
    std::vector<int> pos_max; // stores the end of the allele position, in corresponding index to lines vect
    std::string working_contig;
    bool temp_filled = false;
    bcf1_t * temp_line = bcf_init();
    std::vector<int32_t> temp_dp;
    std::vector<int32_t> temp_ro;
    unsigned short int count = 0;
    if(argc < 4){
        print_help();
    }
    user_input.parse(argc, argv);
    temp_dp.resize((unsigned long) user_input.nsamps);
    temp_ro.resize((unsigned long) user_input.nsamps);
    int32_t * t_dp_dat = temp_dp.data();
    int32_t * t_ro_dat = temp_ro.data();
#if DEBUG
    cerr << "Obtained user inputs, using:\n"
         << "Input:\t" << user_input.vcf_path << "\n"
         << "Reference tag:\t" << user_input.ref_tag << "\n"
         << "Depth tag:\t" << user_input.depth_tag << "\n";
#endif
    while (true){//loop till we can't grab lines
        // clear some things
        std::fill(temp_dp.begin(), temp_dp.end(), INT32_MIN);
        std::fill(temp_ro.begin(), temp_ro.end(), INT32_MIN);
        if (user_input.verbose){
            if(count++ >= user_input.verbose){
                cerr << "At position:\t" << temp_line->pos << "\n";
                count = 0;
            }
        }
        if (temp_filled == false){ // see if we used this line, if not we need a new one
            if(bcf_read(user_input.vcf_file, user_input.hdr, temp_line) != 0){
                if(hts_check_EOF(user_input.vcf_file)){
                    break;
                }
                cerr << "couldn't get a read, but dosn't appere to be end of file....\n";
                break;
            }
            //bcf_unpack(temp_line,BCF_UN_ALL);
        }

#if DEBUG_V
        cerr << "We have a line to work with!\n";
#endif
        if(lines.empty() && pos_max.empty()){
#if DEBUG_V
            cerr << "Our buffer is empty filling in data\n";
#endif
            // it looks like we ether just got rid of our buffer, or are just starting, so put tmp
            // to the first index
            lines.push_back(bcf_init()); // add a spot for our data
            bcf_copy(lines.back(), temp_line); // copy our data in
            pos_max.push_back(temp_line->pos + temp_line->rlen);
            // now to grab the format lines
            bcf_unpack(lines.back(), BCF_UN_FMT);
#if DEBUG_V
            cerr << "Copied in the line \n"
                 << "We have " << user_input.nsamps << " samples!\n";
#endif

            bcf_get_format_int32(user_input.hdr, lines[0], user_input.depth_tag.c_str(),
                                 &user_input.dp_dat, &user_input.nsamps);
            bcf_get_format_int32(user_input.hdr, lines[0], user_input.ref_tag.c_str(),
                                 &user_input.ref_dat, &user_input.nsamps);
            //fix 0 depth, as this should be null for ref and dp
            for (int i = 0 ; i < user_input.depth.size(); i++){
                if(user_input.depth[i] <= 0){
                    user_input.depth[i] = NULL_int32_t;
                    user_input.ref[i] = NULL_int32_t;
                }
            }
            temp_filled = false;
#if DEBUG_V
            cerr <<"At position " << lines[0]->pos << "Obtained tags\nDP-RO:\t";
            for(int i = 0 ; i < user_input.depth.size(); i++){
                cerr << user_input.depth[i] << "-" << user_input.ref[i] << "\t";
            }
            cerr << "\n";
#endif
        } else if (temp_line->pos + temp_line->rlen <= pos_max[0] && temp_line->rid == lines[0]->rid){
#if DEBUG_V
            cerr << "We have data in the range of our buffer, so adding a line\n";
#endif
            // if this is true in a sorted file then we have an ovelap to fix...
            lines.push_back(bcf_init()); // add a spot for our data
            bcf_copy(lines.back(), temp_line); // copy our data in
            pos_max.push_back(temp_line->pos + temp_line->rlen);
            bcf_unpack(lines.back(), BCF_UN_FMT);
#if DEBUG_V
            cerr << "We were able to unpack line, now getting tags!\n";
#endif
            bcf_get_format_int32(user_input.hdr, lines.back(), user_input.depth_tag.c_str(),
                                 &t_dp_dat, &user_input.nsamps);
            bcf_get_format_int32(user_input.hdr, lines.back(), user_input.ref_tag.c_str(),
                                 &t_ro_dat, &user_input.nsamps);
#if DEBUG_V
            cerr << "We were able to obtain tags, now seeing if they should be saved!\n";
#endif
            for (int i = 0; i < user_input.nsamps; ++i) {
                // if the new values are smaller, or not null when storing a null value add them, if not don't
                // ie we are getting non-null min value
                if (temp_dp[i] <= 0) { // if you have a depth of zero, this is wrong, it is null and so is ref...
                    temp_dp[i] = NULL_int32_t;
                    temp_ro[i] = NULL_int32_t;
                } else { // if we jut set things to null don't bother checking
                    if ((temp_ro[i] < user_input.ref[i] && temp_ro[i] != NULL_int32_t)
                            || (user_input.ref[i] == NULL_int32_t )){
                        user_input.ref[i] = temp_ro[i];
                    }
                    if ((temp_dp[i] < user_input.depth[i] && temp_dp[i] != NULL_int32_t)
                            ||(user_input.depth[i] == NULL_int32_t)){
                        user_input.depth[i] = temp_dp[i];
                    }

                }
            }
            temp_filled = false;
#if DEBUG_V
            cerr << "In an overlap block from " << lines[0]->pos << " to " << lines[0]->pos + lines[0]->rlen << std::endl;
            cerr <<"At position " << temp_line->pos << "Updated tags\nDP-RO:\t";
            for(int i = 0 ; i < user_input.depth.size(); i++){
                cerr << user_input.depth[i] << "-" << user_input.ref[i] << "\t";
            }
            cerr << "\n";
#endif
        } else if (!lines.empty() and (
                   temp_line->rid != lines[0]->rid or temp_line->pos + temp_line->rlen > pos_max[0])
                   ){
#if DEBUG_V
            cerr << "We need to output info as\npos:\t" << temp_line->pos << "\nmax_pos:\t" << pos_max[0] << "\n";
#endif
            temp_filled = true;
            output_lines(user_input, lines, pos_max);


        }
    }
    bcf_close(user_input.out_vcf);
    bcf_close(user_input.vcf_file);
    bcf_hdr_destroy(user_input.hdr);
    return 0;
}