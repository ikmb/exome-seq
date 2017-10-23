#!/bin/env ruby
# == NAME
# samplesheet_from_folder.rb
#
# == USAGE
# ./this_script.rb [ -h | --help ]
#[ -i | --infile ] |[ -o | --outfile ] | 
# == DESCRIPTION
# A script to produce a basic sample sheet for exome processing from a folder of fastQ files
#
# == OPTIONS
# -h,--help Show help
# -i,--infile=INFILE input file
# -o,--outfile=OUTFILE : output file

#
# == EXPERT OPTIONS
#
# == AUTHOR
#  Marc Hoeppner, mphoeppner@gmail.com

require 'optparse'
require 'ostruct'


### Define modules and classes here


### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.banner = "Reads Fastq files from a folder and writes a sample sheet to STDOUT"
opts.separator ""
opts.on("-f","--folder", "=FOLDER","Folder to scan") {|argument| options.folder = argument }
opts.on("-h","--help","Display the usage information") {
 puts opts
 exit
}

opts.parse! 

abort "Folder not found (#{options.folder})" unless File.directory?(options.folder)

date = Time.now.strftime("%Y-%m-%d")
center = "ICMB"

fastq_files = Dir["#{options.folder}/*_R*.fastq.gz"]

groups = fastq_files.group_by{|f| f.split("/")[-1].split("_")[0] }

puts "IndivID;SampleID;libraryID;rgID;rgPU;platform;platform_model;Center;Date;R1;R2"

#AS-167616-LR-25610_R1.fastq.gz

groups.each do |group, files|

        left,right = files.sort.collect{|f| File.absolute_path(f)}

        library = group.split("_")[0]
        sample = group.split("-")[0]

        e = `zcat #{left} | head -n1 `

        header = e.to_s.strip.split(":")

        readgroup = "#{header[0]}.#{header[2]}.#{header[3]}.#{library}"

        pgu = header[0].gsub("@", "")

        puts "Indiv_#{sample};Sam_#{sample};#{library};#{readgroup};#{pgu};Illumina;HiSeq3000;#{center};#{date};#{left};#{right}"

end


