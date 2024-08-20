
import os
import glob
from subprocess import Popen
from joblib import Parallel, delayed
from tqdm import tqdm






def estimate_gsize(hap_name):
    adjust_factor = 1
    ps_info = hap_name.split('_')
    # set to be not less than 10k
    gsize = str(max(round(adjust_factor * (int(ps_info[2]) - int(ps_info[1])) / 1000,1), 10)) + "k"
    return gsize

def load_hap_names(filename):
    hap_names = []
    with open(filename, 'r') as file:
        for line in file:
            hap_names.append(line.strip())  # Remove leading/trailing whitespaces     
    return hap_names

def search_fastq_files(hap_names, input_folder):
    fastq_files = []
    for hap_name in hap_names:
        pattern = os.path.join(input_folder, f"/{hap_name}.fastq")
        matching_files = glob.glob(pattern)
        # assert len(matching_files)==1
        fastq_files.extend(matching_files)
    return fastq_files


def reformat_fasta(infile, outfile, hap_name):
    cnt = -1
    with open(infile, 'r') as fin:
        with open(outfile,'w') as fout:
            for line in fin:
                if line[0]=='>':
                    cnt+=1
                    line = ">"+hap_name+'_'+str(cnt)+'\n'
                fout.write(line)
    

def assemble_one_hap_wtdbg2(fq_file, work_dir, out_dir, hap_name,asm_thread, log_file, fail_log_file, data_type, clean ):
    asm_dir = work_dir+"/"+hap_name+"/"
    if not os.path.exists(asm_dir):
        os.system("mkdir -p " + asm_dir)
    else:
        os.system(f"rm -r {asm_dir}; mkdir -p " + asm_dir)

    # asm_len = int((int(hap_name.split('_')[2]) - int(hap_name.split('_')[1]))/1000) + 1

    gsize = estimate_gsize(hap_name)
    # wtdbg2_dir="/data/maiziezhou_lab/Yichen/Softwares/wtdbg2"
    code_dir = os.path.dirname(os.path.realpath(__file__))+"/"
    wtdbg2_dir = code_dir + "/wtdbg2"

    if data_type == 'ONT':
        prefix="ont"
    elif data_type == 'CLR-rs':
        prefix = 'rs'
    else:
        prefix = 'sq'

    cmd = f'''{wtdbg2_dir}/wtdbg2.pl \
        -t {asm_thread} -x {prefix} \
            -g {gsize} -o {asm_dir}/{hap_name} {fq_file}'''
    print(cmd)
    Popen(cmd, shell= True).wait()
    out_file = asm_dir+"/"+hap_name+".cns.fa"

    if os.path.exists(out_file):
        reformat_fasta(out_file, out_dir+"/"+hap_name+".fa" , hap_name)
        cmd = "echo " + hap_name + " >> " + log_file
        Popen(cmd, shell= True).wait()
    else:
        cmd = "echo " + hap_name + " >> " + fail_log_file
        Popen(cmd, shell= True).wait()


    if clean:
        cmd = "rm -r "+ asm_dir
        Popen(cmd, shell= True).wait()

    
def assemble_one_hap_canu(fq_file, work_dir, out_dir, hap_name,asm_thread, log_file, fail_log_file, data_type,clean ):
    asm_dir = work_dir+"/"+hap_name+"/"
    if not os.path.exists(asm_dir):
        os.system("mkdir -p " + asm_dir)
    else:
        os.system(f"rm -r {asm_dir}; mkdir -p " + asm_dir)

    # asm_len = int((int(hap_name.split('_')[2]) - int(hap_name.split('_')[1]))/1000) + 1

    asm_len = estimate_gsize(hap_name)
    # canu_dir="/data/maiziezhou_lab/Yichen/Projects/Long_reads_project/canu-2.1.1/bin"
    code_dir = os.path.dirname(os.path.realpath(__file__))+"/"
    canu_dir = code_dir + "/canu-2.1.1/bin"
    if data_type == 'ONT':
        prefix="-nanopore"
    else:
        prefix = "-pacbio"
    cmd = f'''{canu_dir}/canu -p {hap_name} -d {asm_dir} genomeSize={asm_len} \
        useGrid=false maxThreads={asm_thread} {prefix} {fq_file} '''
    print(cmd)
    Popen(cmd, shell= True).wait()
    out_file = asm_dir+"/"+hap_name+".contigs.fasta"

    if os.path.exists(out_file):
        reformat_fasta(out_file, out_dir+"/"+hap_name+".fa" , hap_name)
        cmd = "echo " + hap_name + " >> " + log_file
        Popen(cmd, shell= True).wait()
    else:
        cmd = "echo " + hap_name + " >> " + fail_log_file
        Popen(cmd, shell= True).wait()


    if clean:
        cmd = "rm -r "+ asm_dir
        Popen(cmd, shell= True).wait()

def assemble_one_hap_miniasm(fq_file, work_dir, out_dir, hap_name,asm_thread, log_file, fail_log_file, data_type,clean ):
    asm_dir = work_dir+"/"+hap_name+"/"
    if not os.path.exists(asm_dir):
        os.system("mkdir -p " + asm_dir)
    else:
        os.system(f"rm -r {asm_dir}; mkdir -p " + asm_dir)


    # miniasm_dir="/data/maiziezhou_lab/CanLuo/Software/miniasm"
    code_dir = os.path.dirname(os.path.realpath(__file__))+"/"
    miniasm_dir = code_dir + "/miniasm"
    if data_type == 'ONT':
        prefix="ont"
    else:
        prefix = 'pb'

    gfa_file = f"{asm_dir}{hap_name}.gfa"
    out_file = asm_dir+"/"+hap_name+".fa"
    cmd = f'''minimap2 -x ava-{prefix} -t{asm_thread} {fq_file} {fq_file} | gzip -1 > {asm_dir}reads.paf.gz;
{miniasm_dir}/miniasm -f {fq_file} {asm_dir}reads.paf.gz > {gfa_file}'''
    print(cmd)
    Popen(cmd, shell= True).wait()

    cmd = '''awk '/^S/{print \">\"$2\"\\n\"$3}' %s | fold > %s '''%(gfa_file, out_file)
    print(cmd)
    Popen(cmd, shell= True).wait()
    
    
    if os.path.exists(out_file):
        reformat_fasta(out_file, out_dir+"/"+hap_name+".fa" , hap_name)
        cmd = "echo " + hap_name + " >> " + log_file
        Popen(cmd, shell= True).wait()
    else:
        cmd = "echo " + hap_name + " >> " + fail_log_file
        Popen(cmd, shell= True).wait()

    if clean:
        cmd = "rm -r "+ asm_dir
        Popen(cmd, shell= True).wait()


def assemble_one_hap_shasta(fq_file, work_dir, out_dir, hap_name,asm_thread, log_file, fail_log_file, data_type, shasta_ont_config, clean ):
    asm_dir = work_dir+"/"+hap_name+"/"
    # shasta does not allow the output dir to be pre-existing
    if os.path.exists(asm_dir):
        os.system(f"rm -r {asm_dir}")


    # shasta_dir="/data/maiziezhou_lab/Yichen/Softwares/Shasta"
    code_dir = os.path.dirname(os.path.realpath(__file__))+"/"
    shasta_dir = code_dir + "/Shasta"

    assert data_type == 'ONT'
    command_file = asm_dir+"command.log"


    cmd = f'''{shasta_dir}/shasta-Linux-0.10.0 \
        --threads {asm_thread} --config {shasta_ont_config} \
            --assemblyDirectory {asm_dir} --input {fq_file}'''
    print(cmd)
    
    Popen(cmd, shell= True).wait()
    with open(command_file,'w') as f:
        f.write(cmd)
        
    out_file = asm_dir+"/Assembly.fasta"

    if os.path.exists(out_file):
        reformat_fasta(out_file, out_dir+"/"+hap_name+".fa" , hap_name)
        cmd = "echo " + hap_name + " >> " + log_file
        Popen(cmd, shell= True).wait()
    else:
        cmd = "echo " + hap_name + " >> " + fail_log_file
        Popen(cmd, shell= True).wait()


    if clean:
        cmd = "rm -r "+ asm_dir
        Popen(cmd, shell= True).wait()

def assemble_one_hap_nextdenovo(fq_file, work_dir, out_dir, hap_name,asm_thread, log_file, fail_log_file, data_type,  clean ):
    asm_dir = os.path.abspath(work_dir+"/"+hap_name+"/")
    nd_out_dir = asm_dir + "/nextdenovo_output/"

    if not os.path.exists(asm_dir):
        os.system(f"mkdir -p {asm_dir}")
    if os.path.exists(nd_out_dir):
        os.system(f"rm -r  {nd_out_dir}")

    if data_type == "ONT":
        dtype_use = "ont"
    else:
        dtype_use = 'clr'
    
    with open(f"{asm_dir}/input.fofn",'w') as f:
        f.write(fq_file+"\n")

    fq_size = os.path.getsize(fq_file)

    # gsize = str(int(fq_size /2/ cov/ 1000) + 1) + "k"

    gsize = estimate_gsize(hap_name)
    code_dir = os.path.dirname(os.path.realpath(__file__))+"/"
    with open(code_dir + "/nextdenovo_config_template.txt",'r') as f:
        s = f.read().\
        replace("<asm_dir>",asm_dir).\
        replace("<ta>",str(asm_thread)).\
        replace("<out_dir>",nd_out_dir).\
        replace("<dtype>", dtype_use).\
        replace("<gsize>",gsize)

    with open(f"{asm_dir}/run.cfg",'w') as f:
        f.write(s)

    # nd_dir="/data/maiziezhou_lab/CanLuo/Software/NextDenovo"
    code_dir = os.path.dirname(os.path.realpath(__file__))+"/"
    nd_dir = code_dir + "/NextDenovo"
    command_file = asm_dir+"command.log"
    cmd = f'''{nd_dir}/nextDenovo {asm_dir}/run.cfg'''
    print(cmd)
    
    Popen(cmd, shell= True).wait()
    with open(command_file,'w') as f:
        f.write(cmd)
        
    out_file = nd_out_dir+"/03.ctg_graph/nd.asm.fasta"

    if os.path.exists(out_file):
        reformat_fasta(out_file, out_dir+"/"+hap_name+".fa" , hap_name)
        cmd = "echo " + hap_name + " >> " + log_file
        Popen(cmd, shell= True).wait()
    else:
        cmd = "echo " + hap_name + " >> " + fail_log_file
        Popen(cmd, shell= True).wait()


    if clean:
        cmd = "rm -r "+ nd_out_dir
        Popen(cmd, shell= True).wait()

def assemble_one_hap_hifiasm(fq_file, work_dir, out_dir, hap_name,asm_thread, log_file, fail_log_file, clean):
    asm_dir = work_dir+"/"+hap_name+"/"
    if not os.path.exists(asm_dir):
        os.system("mkdir -p " + asm_dir)
    else:
        os.system(f"rm -r {asm_dir}; mkdir -p " + asm_dir)

    hifiasm_dir=os.path.dirname(os.path.realpath(__file__))+'/hifiasm-0.14/'
    # hifiasm_dir = "/data/maiziezhou_lab/CanLuo/Software/hifiasm-0.14/"
    # run assembly
    cmd = hifiasm_dir + "/hifiasm -o %s/%s -t %d \
    %s"%(asm_dir,hap_name,asm_thread, fq_file )
    Popen(cmd, shell= True).wait()


    # convert gfa to fa
    cmd = "awk '/^S/{print \">\"$2;print $3}' %s/%s.p_ctg.gfa > %s/%s.p_ctg.fa"%(asm_dir,hap_name,asm_dir,hap_name,)
    Popen(cmd, shell= True).wait()

    # reformat fa
    out_file = "%s/%s.p_ctg.fa"%(asm_dir,hap_name)

    if os.path.exists(out_file):
        reformat_fasta(out_file, out_dir+"/"+hap_name+".fa" , hap_name)
        cmd = "echo " + hap_name + " >> " + log_file
        Popen(cmd, shell= True).wait()
    else:
        cmd = "echo " + hap_name + " >> " + fail_log_file
        Popen(cmd, shell= True).wait()

    if clean:
        cmd = "rm -r "+ asm_dir
        Popen(cmd, shell= True).wait()


def assemble_one_hap_hicanu(fq_file, work_dir, out_dir, hap_name,asm_thread, log_file, fail_log_file, clean):
    code_dir = os.path.dirname(os.path.realpath(__file__))+"/"
    hicanu_dir= code_dir + "/canu-2.1.1/"

    asm_dir = work_dir+"/"+hap_name+"/"
    
    if not os.path.exists(asm_dir):
        os.system("mkdir -p " + asm_dir)
    else:
        os.system(f"rm -r {asm_dir}; mkdir -p " + asm_dir)
    

    gsize = estimate_gsize(hap_name)

    # run assembly
    cmd = hicanu_dir + "/bin/canu -p %s -d %s genomeSize=%s \
        -useGrid=false maxThreads=%d \
    -pacbio-hifi %s"%(hap_name,asm_dir,gsize, asm_thread,fq_file)
    Popen(cmd, shell= True).wait()

    # reformat fa
    out_file = asm_dir+"/"+ hap_name+".contigs.fasta"
    if os.path.exists(out_file):
        reformat_fasta(out_file, out_dir+"/"+hap_name+".fa" , hap_name)
        cmd = "echo " + hap_name + " >> " + log_file
        Popen(cmd, shell= True).wait()
    else:
        cmd = "echo " + hap_name + " >> " + fail_log_file
        Popen(cmd, shell= True).wait()

    if clean:
        cmd = "rm -r "+ asm_dir
        Popen(cmd, shell= True).wait()

def assemble_one_hap_flye(fq_file, work_dir, out_dir, hap_name,asm_thread, log_file, fail_log_file, data_type, clean):
    asm_dir = work_dir+"/"+hap_name+"/"
    if not os.path.exists(asm_dir):
        os.system("mkdir -p " + asm_dir)
    else:
        os.system(f"rm -r {asm_dir}; mkdir -p " + asm_dir)

    if data_type =='ONT':
        cmd = "flye --nano-raw %s --out-dir %s --threads %d"%(fq_file,asm_dir,asm_thread)
    elif 'CLR' in data_type:
        cmd = "flye --pacbio-raw %s --out-dir %s --threads %d"%(fq_file,asm_dir,asm_thread)
    else:
        print("flye only support ONT and CLR")
        exit()
    
    Popen(cmd, shell = True).wait()

    # reformat fa
    out_file = asm_dir+"/assembly.fasta"
    if os.path.exists(out_file):
        reformat_fasta(out_file, out_dir+"/"+hap_name+".fa" , hap_name)
        cmd = "echo " + hap_name + " >> " + log_file
        Popen(cmd, shell= True).wait()
    else:
        cmd = "echo " + hap_name + " >> " + fail_log_file
        Popen(cmd, shell= True).wait()

    if clean:
        cmd = "rm -r "+ asm_dir
        Popen(cmd, shell= True).wait()


def remove_duplicate_one_fastq(input_path,output_path, log_file):
    fw = open(output_path,'w')
    name_dc = {}
    write_or_not = []
    cnt = 0
    for line in open(input_path):
            if cnt%4==0:
                    name = line[1:-1]
                    if name not in name_dc:
                            name_dc[name] = 1
                            write_or_not.append(1)
                    else:
                            write_or_not.append(0)
            read_num = cnt//4
            decision = write_or_not[read_num]
            if decision==1:
                    fw.write(line)
            cnt+=1
    fw.close()

    os.system( f"echo {input_path} >> {log_file}")

def remove_duplicate(in_list, n_thread, ):

    in_dir_list = list(set([os.path.dirname(in_file) for in_file in in_list]))
    if len(in_dir_list)==0:
        return
    assert len(set(in_dir_list)) == 1
    in_dir  = in_dir_list[0]
    out_dir = in_dir + "_rmdp"
    if not os.path.exists(out_dir):
        os.system("mkdir -p " + out_dir)

    log_file = f"{in_dir}_rmdp.log"
    os.system("touch " + log_file)

    with open(log_file,'r') as f:
        done_files = f.read().split('\n')[:-1]

    pending_files = list(set(in_list) - set(done_files))

    print("In toal %d fastqs"% len(in_list))
    print("%d fastq(s) have already been processed"% len(done_files))
    print("%d fastq(s) still need to be processed"% len(pending_files))


    out_files = [ out_dir + "/" +os.path.basename(fq)   for fq in pending_files ]

    sequences = Parallel(n_jobs=n_thread)(delayed(remove_duplicate_one_fastq)(pending_files[i], out_files[i],log_file ) for i in tqdm(range(len(pending_files ))))




def assemble_one_hap(assembler,fq_file, work_dir, out_dir, hap_name,asm_thread, log_file, fail_log_file, data_type, shasta_ont_config, clean ):

    if assembler == 'wtdbg2':
        assemble_one_hap_wtdbg2(fq_file, work_dir, out_dir, hap_name,asm_thread, log_file, fail_log_file, data_type, clean )
    elif assembler == 'canu':
        assemble_one_hap_canu(fq_file, work_dir, out_dir, hap_name,asm_thread, log_file, fail_log_file, data_type, clean )
    elif assembler == 'miniasm':
        assemble_one_hap_miniasm(fq_file, work_dir, out_dir, hap_name,asm_thread, log_file, fail_log_file, data_type, clean )
    elif assembler == 'shasta':
        assemble_one_hap_shasta(fq_file, work_dir, out_dir, hap_name,asm_thread, log_file, fail_log_file, data_type, shasta_ont_config, clean )
    elif assembler == 'nextdenovo':
        assemble_one_hap_nextdenovo(fq_file, work_dir, out_dir, hap_name,asm_thread, log_file, fail_log_file, data_type,  clean )
    elif assembler == "hifiasm": 
        assemble_one_hap_hifiasm(fq_file, work_dir, out_dir, hap_name,asm_thread, log_file, fail_log_file, clean)
    elif assembler == 'hicanu':
        assemble_one_hap_hicanu(fq_file, work_dir, out_dir, hap_name,asm_thread, log_file, fail_log_file, clean)
    elif assembler == 'flye':
        assemble_one_hap_flye(fq_file, work_dir, out_dir, hap_name,asm_thread, log_file, fail_log_file, data_type, clean)
    else:
        print("only support wtdbg2, canu, miniasm, shasta, nextdenovo, hifiasm, hicanu, flye")
        exit()



# hap_file = args.hap_file
# fastq_dir = args.fastq_dir
# output_dir = args.output_dir
# assembler = args.assembler
# n_thread = args.n_thread
# asm_thread = args.asm_thread
# clean = args.clean
# data_type = args.data_type
# shasta_ont_config = args.shasta_ont_config
# pacbio_subtype = args.pacbio_subtype
# prefix=args.prefix

def run_assembly_one_folder(hap_file, fastq_dir,  use_dir,assembler, 
                 n_thread, asm_thread, clean, data_type, shasta_ont_config,
                  pacbio_subtype, prefix ):
    
    if (assembler == 'shasta') and (shasta_ont_config is None):
        print("You must specify shasta_ont_config when you are using shasta!")
        exit()
    if (assembler == 'shasta') and (data_type != 'ONT'):
        print("Shasta only support ONT data!")
        exit()

    if (data_type == 'Hifi') and (assembler not in ['hifiasm','hicanu']):
        print("Only hifiasm and hicanu support Hifi data assembly!")
        exit()


    global code_dir
    code_dir = os.path.dirname(os.path.realpath(__file__))+'/'

    if pacbio_subtype is not None:
        data_type = data_type + '-'  + pacbio_subtype



    if hap_file is not None:
        # ------- match fastq files
        hap_names = load_hap_names(hap_file)
        print("Num of haplotype names: ",len(hap_names))
        print("example loaded haplotype names:")
        print(hap_names[0] )

        found_files = search_fastq_files(hap_names, fastq_dir)
        print("Num of found fq files: ",len(found_files))
        print("example found FASTQ files:")
        print( found_files[0])

        assert len(found_files) == len(hap_names)
    else:
        pattern = os.path.join(fastq_dir, "*.fastq")
        print(pattern)
        found_files = glob.glob(pattern)
        hap_names = [fq.split('/')[-1].split('.')[0] for fq in found_files]
        assert len(hap_names) == len(set(hap_names))
        print("Num of fastqs: ", len(hap_names))

    
    if len(found_files) == 0:
        print(f"No fastq files are found in {fastq_dir}")
        return
    #---------- remove redundancy
    fq_dir_name = os.path.basename(os.path.normpath(fastq_dir))
    if data_type!= 'Hifi':

        remove_duplicate(found_files, n_thread)
        # update fq list
        found_files = [os.path.dirname(fq)+'_rmdp/'+os.path.basename(fq)   for fq in found_files ]

        output_dir = use_dir + "/Assembly_Output_" + fq_dir_name+'_rmdp'
    else:
        output_dir = use_dir + "/Assembly_Output_" + fq_dir_name
    # ----------- assembly
    

    work_dir = output_dir+"/assembly_temp_files/"
    contig_dir = output_dir+"/fasta_by_hap/"
    final_contig_dir = output_dir + "/final_contigs/"
    if not os.path.exists(work_dir):
        os.system("mkdir -p " + work_dir)

    if not os.path.exists(contig_dir):
        os.system("mkdir -p " + contig_dir)



    log_file = output_dir+"/log.txt"
    os.system("touch " + log_file)

    fail_log_file = output_dir+"/fail_log.txt"
    os.system("touch " + fail_log_file)

    # global cov
    # cov = estimate_coverage(fastq_dir)
    # check previously done asm
    finished_haps = set(load_hap_names(log_file) + load_hap_names(fail_log_file))

    rest_haps = []
    rest_fqs = []
    for i in range(len(hap_names)):
        if hap_names[i] not in finished_haps:
            rest_haps.append(hap_names[i])
            rest_fqs.append(found_files[i])
    print("Num of remaining tasks: ", len(rest_haps))



    sequences = Parallel(n_jobs=n_thread)(delayed(assemble_one_hap)(
        assembler,
        rest_fqs[i], work_dir, 
        contig_dir, rest_haps[i],
        asm_thread, log_file, 
        fail_log_file, data_type, 
        shasta_ont_config, clean ) 
                                        for i in range(len(rest_haps)))


    # collect all
    if not os.path.exists(final_contig_dir):
        os.system("mkdir -p " + final_contig_dir)

    cmd = f"cat {contig_dir}/*.fa > {final_contig_dir }/{prefix}_final_contigs.fa"
    Popen(cmd, shell= True).wait()

def run_assembly(hap_file, fastq_dirs, output_dir, assemblers, 
                 n_thread, asm_thread, clean, data_type, shasta_ont_config,
                  pacbio_subtype, prefix ):
    
    if len(assemblers) == 1:
        assemblers = assemblers * len(fastq_dirs)
    
    assert len(fastq_dirs) == len(assemblers)
    fasta_list = []
    print(fastq_dirs)
    for i in range(len(fastq_dirs)):
        fastq_dir = fastq_dirs[i]
        assembler = assemblers[i]
        print(fastq_dir)
        print(f"\n**********Process {fastq_dir} using {assembler}...\n")
        run_assembly_one_folder(hap_file, fastq_dir,  output_dir, assembler, 
                 n_thread, asm_thread, clean, data_type, shasta_ont_config,
                  pacbio_subtype, prefix )
        if data_type == 'Hifi':
            fasta = os.path.dirname(os.path.normpath(fastq_dir)) + "/Assembly_Output_" + os.path.basename(os.path.normpath(fastq_dir))+f"/final_contigs/{prefix}_final_contigs.fa"
        else:
            fasta = os.path.dirname(os.path.normpath(fastq_dir)) + "/Assembly_Output_" + os.path.basename(os.path.normpath(fastq_dir))+f"_rmdp/final_contigs/{prefix}_final_contigs.fa"
        fasta_list.append(fasta)
        
    ##### concat all fasta
    final_contig_dir = output_dir+ "/final_contigs/"
    if not os.path.exists(final_contig_dir):
        os.system("mkdir -p " + final_contig_dir)

    cmd = f"cat {' '.join(fasta_list)} > {final_contig_dir}/{prefix}_final_contigs.fa"
    Popen(cmd, shell = True).wait()
        
if __name__ == "__main__":
    import argparse
    from argparse import ArgumentParser
    parser = ArgumentParser(description="",
        usage='use "python3 %(prog)s --help" for more information',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--hap_file','-haps')
    parser.add_argument('--fastq_dirs','-fqds',
                        help='the folders that includes FASTQ files.', nargs = '+')
    parser.add_argument('--output_dir','-o')
    parser.add_argument('--assemblers','-asms', choices = ['wtdbg2','canu','miniasm','shasta','nextdenovo','hifiasm','hicanu','flye'], nargs= '+', 
                        help = "the assemblers used for fastq_dirs; order corresponds to the order of fastq_dirs")
    parser.add_argument('--data_type','-d', choices = ['CLR','ONT','Hifi'])
    parser.add_argument('--pacbio_subtype','-pb', choices = ['rs','sq'], help = "must provide when using wtdbg2 on CLR data")
    parser.add_argument('--shasta_ont_config','-shacon', choices = ['Nanopore-OldGuppy-Sep2020'])
    parser.add_argument('--n_thread','-t', type = int, default = 7 )
    parser.add_argument('--asm_thread','-ta', type = int, default = 7 )
    parser.add_argument('--clean','-cl', action='store_true')
    parser.add_argument('--prefix','-px', help = "file prefix in the output folder", default = "Sample")
    args = parser.parse_args()
    global data_type, prefix
    hap_file = args.hap_file
    fastq_dirs = args.fastq_dirs
    output_dir = args.output_dir
    assemblers = args.assemblers
    n_thread = args.n_thread
    asm_thread = args.asm_thread
    clean = args.clean
    data_type = args.data_type
    shasta_ont_config = args.shasta_ont_config
    pacbio_subtype = args.pacbio_subtype
    prefix=args.prefix
    run_assembly(hap_file, fastq_dirs, output_dir, assemblers, 
                 n_thread, asm_thread, clean, data_type, shasta_ont_config,
                  pacbio_subtype, prefix )
