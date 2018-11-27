import pathlib
import os
import pandas
import jinja2
import sh
import logging
import filecmp
import datetime
import numpy
import itertools
import subprocess
from cleo import Command, argument, option
from Bio import SeqIO, Phylo



class RunCoreugate(Command):

    '''
    A class to set up file structure for running Coreugate assembly/chewBBACA pipeline
    '''
    # name of command to run
    name = 'run'

    arguments = [

    ]
    # user defined or modifable options
    options = [
        option('input_file', 'f', description= 'Input file - three columns <isolate_id> <path/to/R1> <path/to/R2>', default='', value_required=True),
        option('schema', 'd', description='Path to species schema/allele db', default='', value_required= True),
        option('workdir', 'w', description = 'Working directory, defaults to CWD' ,default=f"{pathlib.Path.cwd().absolute()}", value_required=True),
        option('template_dir', 'r', description='Path to Snakefile and config.yaml templates.', default= f"{pathlib.Path(__file__).parent.parent / 'templates'}", value_required=True),
        option('threads', 't', description='Number of threads for chewBBACA', default=36, value_required=True),
        option('id', 'i', description= 'Name of the job', default = '' ,value_required=True),
        option('singularity', 's', description='Use singularity containers for assembly and chewBBACA rather than local install', default='N', value_required=True),
        option('min_contig_size', 'c', description='Minumum contig size required for QC', default=500, value_required=True),
        option('min_contigs', 'm', description='Minumum number of contigs required for QC', default=0, value_required=True),
        option('assembler', 'a', description='Assembler to be used (options are: shovill-spades, shovil-skesa, shovill-velvet, skesa, spades)', default= 'shovill-spades', value_required=True),
        option('prodigal_training', 'p', description='Prodigal file to be used in allele calling. See https://github.com/B-UMMI/chewBBACA/tree/master/CHEWBBACA/prodigal_training_files for options', default= '', value_required=True),
        
    ]

    def handle(self):
        # check singularity and versions
        singularity = f"{self.option('singularity')}".upper()
        assembler = self.option('assembler')
        if singularity == 'Y':
            self.log_messages('info', 'You have selected to run Coreugate with Singularity containers')
        elif singularity == 'N':
            self.check_chewbbaca()
            self.check_assemblers()
                   
        # get the datetime for labels where needed
        now = datetime.datetime.today().strftime("%d_%m_%y_%H")
        day = datetime.datetime.today().strftime("%d_%m_%y")
        # check if templates and working dir are present
        template_dir = pathlib.Path(self.option('template_dir'))
        self.path_exists(template_dir)
        workdir = pathlib.Path(self.option('workdir'))
        self.path_exists(workdir)
        
        # check input file       
        if self.option('input_file') == '':
            ifile = self.ask('No input file added. Please provide path to input file.')
        else:
            ifile = self.option('input_file')
        input_file = pathlib.Path(ifile)
        self.path_exists(input_file)
        
        
        # check if id has been given
        job_id = self.option('id')
        if self.id_exists(job_id) == False:
            job_id = self.ask('Job ID can not be empty. Please provide a unique job ID.')
        # threads for chewbacca
        cpu = self.option('threads')
        # set up schema and wroking directory
        schema = pathlib.Path(self.option('schema'))
        # assemble is Y or N depending the data input types, job_directory is a path object and schema path a path object to linked schem path
        assemble, job_directory, schema_path = self.setup_working_directory(input_path = input_file, workdir = workdir, day= day, job_id= job_id, schema = schema, cpu = cpu)
        
        min_contig_size = self.option('min_contig_size')
        min_contigs = self.option('min_contigs')
        # set ptf option for chewbbaca - if empty then warn user
        ptf = f"{self.option('prodigal_training')}"
        if ptf == '':
            self.log_messages('info', 'You are running chewBBACA without a training file. Please note this is against developers recomendations.')

        # setup Snakefile and config file
        self.write_workflow(workdir=workdir, template_dir=template_dir,assemble=assemble,assembler=assembler, min_contig_size=min_contig_size, min_contigs=min_contigs, cpu=cpu, schema_path=schema_path, job_directory=job_directory, ptf = ptf)
        # run the pipeline
        if self.run_workflow(job_directory = job_directory, singularity = singularity) == True:
            self.finish_workflow()

        
    def check_version(self, software):
        '''
        check version of software installed, if chewbbaca add string
        :software: name of software (str)
        '''
        sft = subprocess.run([software, '--version'], stdout=subprocess.PIPE)
        if software == 'chewBBACA.py':
            version = f"chewBBACA v{sft.stdout.decode().strip()}"
        else:
            version = f"{sft.stdout.decode().strip()}"
        return version
        
    def check_chewbbaca(self):
        '''
        check chewbbaca version -> needs to be 2.0.16
        '''
        try:
            chewie_version = self.check_version('chewBBACA.py')
            if chewie_version == f"chewBBACA v2.0.16":
                self.log_messages('info', f"chewBBACA version {chewie_version} has been found.")
            else:
                self.log_messages('warning', f"{chewie_version} has been found. This is not compatible with COREugate, please install chewBBACA version 2.0.16 before proceeding (https://github.com/B-UMMI/chewBBACA). Or use `-s Y`. Exiting....")
                raise SystemExit
        except FileNotFoundError:
            self.log_messages('warning', f"chewBBACA is not installed. please install chewBBACA <= version 2.0.16 before proceeding (https://github.com/B-UMMI/chewBBACA). Or use `-s Y`. Exiting....")
            raise SystemExit
    
    def check_velvet(self):
        '''
        check velveth and velvetg installation
        '''
        vh = subprocess.run(['velveth', '--version'], stdout=subprocess.PIPE)
        vg = subprocess.run(['velvetg', '--version'], stdout=subprocess.PIPE)
        if vh.returncode == 0 and vg.returncode == 1:
            return True
        else:
            return False

        
    def check_assemblers(self):
        '''
        check assemblers
        '''
        sftw_dict = {'skesa': 'SKESA v.2.3.0', 'shovill' : 'shovill 1.0.4', 'spades.py' : 'SPAdes v3.13.0'}
        nf = []
        for s in sftw_dict:
            asmb_version = self.check_version(s)
            if asmb_version == sftw_dict[s]:
                self.log_messages('info', f"{asmb_version} has been found.")
            else:
                nf.append(sftw_dict[s])
        if self.check_velvet() == False:
            nf.append('Velvet')
        if len(nf) > 0:
            for n in nf:
                self.log_messages('warning', f"Some software is not compatible with COREugate. Please install {n} or use `-s Y`")
            raise SystemExit

    
    def log_messages(self, type, message):
        '''
        so as to not repeat myself too often...
        '''
        if type == 'warning':
            logging.warning(message)
            print(f"WARNING: {message}")
        if type == 'info':
            logging.info(message)
            print(f"{message}")


    def path_exists(self, path):
        '''
        input: a pathlib PosixPath
        output: if path exists output path is found else raise FileNotFoundError
        '''

        if not path.exists():
            self.log_messages('warning', f"The {path.name} does not exist.")
            raise FileNotFoundError(f"{path.name}")
        else:
            return True

    def id_exists(self, name):
        '''
        check if the name is an empty string ID can not be empty
       
        '''
        if isinstance(name, str):
            if len(name) == 0:
                return False
            else:
                return True
        else:
            return False
    
    def link(self, workdir, path):
        '''
        link is for directory NOT for single files!
        '''
        target = workdir / path.parts[-1]
        if not target.exists():
            self.log_messages('info', f"Linking {path.name} to {workdir.name}")
            target.symlink_to(path)
        
        return True

    def link_schema(self, workdir, path):
        '''
        create a link to the external schema in the working dir
        '''

        found = False
        while not found: 
            
            if self.path_exists(path=path):
                
                found = self.link(workdir = workdir, path = path)
                
            else:
                path = self.ask(f"Path to schema does not exist. Please enter a valid path: ")
                path = pathlib.Path(path)
        
        return  path.name

    
    def prep_external_schema(self, path, cpu):
        '''
        check if the schema has been prepped
        :path:
            path is path to schema
        :cpu:
            cpu for running PrepExternalSchema (threads option)
        '''
        short_path= path / 'short'
        # if no short directory run PrepExternalSchema
        if not os.path.exists(short_path):
            print('Preparing external schema')
            subprocess.call(['chewBBACA.py', 'PrepExternalSchema', '-i', f"{path}", '--cpu', f"{cpu}", '-v' ])
        else:
            self.log_messages('info','Schema is already in the correct format. Congratulations.')
        # if short/fasta is present and empty then remove it to prevent problems with chewbbaca
        fasta = path.glob("*.fasta")
        self.log_messages('info', 'Checking for abberent fasta files. Incorrectly formatted fasta files and empty fasta files will be removed.')
        for f in fasta:
            if os.path.exists(short_path / f / '_short.fasta') and os.path.getsize(short_path / f / '_short.fasta') == 0:
                subprocess.call(['rm', short_path / f / '*'])
                subprocess.call(['rm', path / f / '*'])
        # if temp exists remove it
        if os.path.exists(path / 'temp'):
            subprocess.call(['rm', '-r',  str(path/ 'temp')])
        
    
    def check_input_exists(self, df, data_type):
        '''
        check that the correct paths have been given
        if reads/assemblies not present path_exists will cause a FileNotFound error and warn user
        :df: a dataframe of user inpt file
        :data_type:'reads' or 'assemblies' decide whether there should be 2 or 3 
        '''
        
        for i in df.itertuples():
            r1 = i[2]
            self.path_exists(pathlib.Path(r1))
            if data_type == 'reads':
                r2 = i[3]
                self.path_exists(pathlib.Path(r2))
        return True
    
    
    def set_isolate_log(self, isolate_df, workdir, day, data_type, job_id):
        '''
        add the isolates to a log file, if there is already a logile of the given datatype then open it check for changes and save
        input:
            :isolate_df: dataframe of the isolates to add 
            :workdir: path to working directory (where log will/is be)
            :day: day of the run for log (datetime.datetime  as string)
            :logfile: path to logfile which will have the datatype in the name - if user wants to add assemblies next time - although this is perhaps not advisable
        '''        
        
        logfile = workdir / job_id / f"job_{data_type}.log"
        # check that input exists - if not inform user and quit.
        self.check_input_exists(isolate_df, data_type)
        # if assemblies df has two columns
        if data_type == 'assemblies':
            lf = pandas.DataFrame({'Isolate': [i for i in list(isolate_df.iloc[ : , 0]) if '#' not in i ], 'Assemblies' : [i for i in list(isolate_df.iloc[ : , 1])],'Status': f"INCLUDED", 'Date': day})
        # else has 3
        else:
            lf = pandas.DataFrame({'Isolate': [i for i in list(isolate_df.iloc[ : , 0]) if '#' not in i ], 'R1' : [i for i in list(isolate_df.iloc[ : , 1])], 'R2': [i for i in list(isolate_df.iloc[ : , 2])], 'Status': f"INCLUDED", 'Date': day})
        # if there is already a logfile for that data type open it
        if logfile.exists():            
            old_df = pandas.read_csv(logfile, sep = '\t', index_col = False)
            # it there are differences append the new isolates and remove duplications
            if old_df.equals(lf) == False:
                lf = old_df.append(lf)
                lf = lf.drop_duplicates()
            else:
                lf = old_df
        # save the file 
        lf.to_csv(logfile, sep = '\t', index = False)    
        
        return True
     

    def link_files(self, target, source):
        '''
        for linking reads and assemblie to the data directory
        '''
        try:
            target.symlink_to(source)
        except FileExistsError:
            self.log_messages('warning', f"{target} already exists.")

        # else:
        #     self.log_messages('info', f"{target} already exists.")

    def link_reads(self, read_name, data_dir_name, row):
        '''
        take read_name, data_dir_name (assemblies or READS) and path to source (row) to construct a target and source path
        ''' 
        # from the input file df - already checked for existence
        source = pathlib.Path(row)
        # where the file will be put
        target = data_dir_name / read_name
        # make sure the data directory exists - otherwise this will fail
        if not data_dir_name.exists():
            data_dir_name.mkdir()
        # link
        if not target.exists():
            self.link_files(target = target , source = source)
        return True

    def make_dir(self, data_dir):
        '''
        make a data directory (either assemblies or READS)
        '''
        if not data_dir.exists():
            data_dir.mkdir()
    
    def reads_or_assembly(self,df):
        '''
        check whether the data is likely to be assemblies or reads based on input dimensions.
        '''
        if df.shape[1] ==3:
            return 'reads'
        elif df.shape[1] == 2:
            return 'assemblies'
        else:
            return False  
    
    
    def set_input_file(self, path, workdir, day, job_id):
        '''
        setup the input file - checks the type of input data, and checks if data is available
        '''
        assemble = 'N' # Default N, if there are 3 cols and reads_or_assembly returns reads, this will change to Y
        # check if using reads or assemblies
        df = pandas.read_csv(path, sep = None,engine = 'python', header = None)
        data_type = self.reads_or_assembly(df)

        if data_type == False:
            self.log_messages('info', 'The input file is not in the correct format and can not be read. Please check and try again.')
            raise SystemExit
        else:
            self.log_messages('info',f"You have elected to use {data_type} as an input source.")
            if data_type == 'reads':
                assemble = 'Y'
        # set up log file will check and confirm that all input files exists
        if self.set_isolate_log(isolate_df= df, workdir = workdir, day = day, data_type = data_type, job_id=job_id):
            df.to_csv(workdir / 'isolates.tab', sep = '\t', header = False, index = False)
            self.log_messages('info', 'Successfully saved `isolates.tab`')
        
        return assemble

    def setup_working_directory(self, input_path, workdir, day, job_id, schema, cpu):
        '''
        ensure all required files are linked to the working directory, job_directory and data directory
        :input_path: input_file file path
        :workdir: working directory
        :day: day for log
        :job_id: unique identifer will be used to make the job_directory
        :schema: path to schema
        :cpu: cpu for running PrepExternalSchema if needed
        '''
        # make the job directory
        job_directory = workdir / f"{job_id}"
        if not job_directory.exists():
            job_directory.mkdir()
        
        # link schema to working directory
        schema_path = pathlib.Path(self.link_schema(workdir = job_directory, path = schema))
        self.prep_external_schema(path = job_directory / schema_path, cpu = cpu)

        # if assemble = Y then the input directory is reads, otherwise assemblies
        # set_input_file will generate the isolates.tab file for running and also the
        assemble = self.set_input_file(path = input_path, workdir=workdir, day=day, job_id= job_id)
        
        # get input file
        isolates_tab = pandas.read_csv(workdir / 'isolates.tab', sep = '\t', index_col = False,header = None)
        # set up data directory
        if assemble == 'Y':
            data_dir = job_directory / 'READS'
        else:
            data_dir = job_directory / 'assemblies'
        # make the data_dir for linking into
        self.make_dir(data_dir = data_dir)
        # link input data to input directory
        self.log_messages('info', f"Linking source data to working directory")
        for row in isolates_tab.itertuples():
            # for assemblies
            name = row[1]
            if assemble == 'N':
                source = pathlib.Path(row[2]) # source has already been confirmed as exisiting no need to re-check
                target = data_dir / f"{name}.fa" # name of isolate .fa
                self.link_files(target = target, source=source)
            else:
                data_dir_name = data_dir / name
                

                self.link_reads(data_dir_name=data_dir_name, read_name= 'R1.fq.gz', row = row[2])
                self.link_reads(data_dir_name=data_dir_name, read_name= 'R2.fq.gz', row = row[3])

        return assemble, job_directory, schema_path

    
    def assembly_string(self):
        '''
        the string for insertion into the Snakefile if assembly is needed
        '''

        assembly_string = f"""
rule assembly:
    input:
        r1 = \"READS/{{sample}}/R1.fq.gz\",
        r2 = \"READS/{{sample}}/R2.fq.gz\"
    output:
        \"assemblies/{{sample}}.fa\"
    singularity:
        \"shub://phgenomics-singularity/multi_assembler_singularity@latest\"

    shell:
        \"""
        if [ \"{{assembler}}\" == "shovill-spades" ]; then
            shovill --outdir temp --R1 {{input.r1}} --R2 {{input.r2}} --force --minlen 200
            cp temp/contigs.fa {{output}}
            rm -r temp
        elif [ \"{{assembler}}\" == \"shovill-skesa\" ]; then
            shovill --outdir temp --R1 {{input.r1}} --R2 {{input.r2}} --force --minlen 200 --assembler skesa --opts \"--vector_percent 1\"
            cp temp/contigs.fa {{output}}
            rm -r temp
        elif [ \"{{assembler}}\" == \"shovill-velvet\" ]; then
            shovill --outdir temp --R1 {{input.r1}} --R2 {{input.r2}} --force --minlen 200 --assembler velvet 
            cp temp/contigs.fa {{output}}
            rm -r temp
        elif [ \"{{assembler}}\" == \"skesa\" ]; then
            skesa --fastq {{input.r1}},{{input.r2}} --vector_percent 1 --use_paired_ends --cores 4 > {{output}}
        elif [ \"{{assembler}}\" == \"spades\" ]; then
            spades.py -o temp/ -1 {{input.r1}} -2 {{input.r2}}
            cp temp/contigs.fasta {{output}}
            rm -r temp
        fi
        \"""
"""
        return assembly_string

    def write_workflow(self, workdir, template_dir, assemble, assembler, min_contig_size, min_contigs,cpu,schema_path, job_directory, ptf):
        '''
        write the Snakefile and config.yaml
        '''
        if assemble == 'Y':
            assemble_string = self.assembly_string()
        else:
            assemble_string = ''
        self.log_messages('info',f"Setting up specific workflow")        
        # read the config file which is written with jinja2 placeholders (like django template language)
        config_template = jinja2.Template(pathlib.Path(template_dir, 'config.yaml').read_text())
        config = job_directory / 'config.yaml'
        
        config.write_text(config_template.render(assembler = assembler, assemble = assemble, min_contig_size = min_contig_size, min_contigs = min_contigs,cpu=cpu, schemPath = schema_path,workdir = job_directory))
        
        self.log_messages('info',f"Config file successfully created")

        snk_template = jinja2.Template(pathlib.Path(template_dir, 'Snakefile').read_text())
        snk = job_directory / 'Snakefile'
        snk.write_text(snk_template.render(ptf = f"--ptf {ptf}", assemble_string = assemble_string))
        
        self.log_messages('info',f"Snakefile successfully created")

        rfile_template = jinja2.Template(pathlib.Path(template_dir, 'distances.R').read_text())
        rfile = job_directory / 'distances.R'
        rfile.write_text(rfile_template.render(script_dir = f"{template_dir}"))




    def run_workflow(self, job_directory, singularity):
        '''
        run Coreugate
        set the current directory to working dir for correct running of pipeline if singularity = Y then run with singularity
        if the pipeline wroks, return True else False
        '''

        os.chdir(job_directory)
        
        if singularity == 'Y':
            cmd = f"snakemake -s Snakefile --use-singularity --singularity-args '--bind /home'"
        else:
            cmd = f"snakemake -s Snakefile"
        self.log_messages('info', f"Running {cmd}")
        wkf = subprocess.run(cmd, shell = True)
        if wkf.returncode == 0:
            return True
        else:
            return False
        
    def finish_workflow(self):
        '''
        final message at completion of workflow. If workflow goes to completion print 'thanks for coming message'
        '''
        self.log_messages('info', f"COREugate has finished.")
        self.log_messages('info', f"Have a nice day. Come back soon.") 
        self.log_messages('info',f"{60 * '='}")