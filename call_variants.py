from subprocess import check_output, STDOUT, CalledProcessError

def call_command(command, command_name, cwd=None):
    print("About to invoke {} with command {}.".format(command_name, command)) 
    try:
        my_output = check_output(command, shell=True, cwd=cwd, stderr=STDOUT)
        print(my_output)
    except CalledProcessError as e:
        print("Program output is: " + e.output.decode("utf-8") )
        print("{} execution failed {}.".format(command_name, e.returncode))
        raise

def get_output_location(command_name, sample_name):
    return output_location + command_name + "/" + current_sample_name + "." + command_name + ".vcf"

reference_location = "~/data/human_g1k_v37.20.21.fasta"
region = "20:9999900-10001100"

output_location = "~/data/snp_calling_tests/"

input_location = "~/data/snp_calling_tests/inputs/"
sample_names = ["CEUTrio.HiSeq.WGS.b37.NA12878.20.10000117","CEUTrio.HiSeq.WGS.b37.NA12878.20.10000758"]
sample_suffices = [["01","02","05","10","20","30"], ["01","02","05","10","20"]]

callers = ["freebayes", "samtools", "platypus", "gatk4"]
caller_commands = {"freebayes" : "freebayes -f {} -r {} --report-genotype-likelihood-max {} > {}",
                   "samtools": "samtools mpileup -ugf  {} -r {} {} | bcftools call -vmO v -o {}",
                   "platypus": "python ~/tools/platypus/bin/Platypus.py callVariants --refFile={} --regions={} --bamFiles={} --output={}",
                   "gatk4": "~/tools/gatk4/gatk-protected/gatk-launch HaplotypeCaller -R {} -L {} -I {} -O {} -pairHMM AVX_LOGLESS_CACHING"}



for i in range(len(sample_names)):
    for sample_suffix in sample_suffices[i]:
        current_sample_name = sample_names[i] + "." + str(sample_suffix)
        input_filename = input_location + current_sample_name + ".bam"

        for caller_name in callers:
            output_filename = get_output_location(caller_name, current_sample_name)
            full_command = caller_commands[caller_name].format(reference_location, region, input_filename, output_filename)
            call_command(full_command, caller_name)
            call_command("bgzip " + output_filename, "bgzip")
            call_command("tabix " + output_filename + ".gz", "tabix")
