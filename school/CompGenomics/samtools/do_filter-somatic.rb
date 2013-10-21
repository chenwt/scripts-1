#!/usr/bin/env ruby

require 'getoptlong'

### 

## filter vcf from samtools
# somatic mutations:
# number of alternative alleles from + and - strand: both >= minStrandReads
# phred-scale likelihood: in cancer:  0/0 > minPhredT, 0/1 = 0 or 1/1 = 0; in normal: 0/0 = 0, 0/1 and 1/1 > minPhredN
# PV4 > minPvalue
# mq0 ratio < maxMQ0Ratio



def main 
  settings = {}
  settings["--dp4"]=2
  settings["--pv41"]=0.000001
  settings["--pv42"]=0.000001
  settings["--pv43"]=0.000001
  settings["--pv44"]=0.000001
  settings["--clr"]=20
  settings["--normal"]=0  ## default: tumor is first, normal is second. 
  settings["--minPL"] = 40
  optHash = getopt()
  vcf = optHash["--vcf"]
  
  settings.keys.sort.each do |s|
    if optHash.key?(s)
      settings[s] = optHash[s].to_f
    end
  end
  
  nsample=countSamples(vcf)
  
  filterVCF(vcf,settings,nsample)  # gt: gene -> pos -> sample -> genotype, 

end

def filterVCF(vcf, settings, nsample)
  o = File.new(vcf + ".somaticfiltered.vcf", 'w')
  e = File.new(vcf + ".dropped.vcf", 'w')
  firstline = 1
  fields = {}

  File.new(vcf, 'r').each do |line|
    if line.match("^#")
      o.puts line
      e.puts line
    else
      cols=line.chomp.split(/\t/)
      qual, pass, info = cols[5].to_f, cols[6], cols[7].split(';')
      #10      61852   .       T       C       999     .       DP=76;AF1=0.5;AC1=2;DP4=9,17,11,15;MQ=40;FQ=999;PV4=0.78,0.17,0.42,0.043        GT:PL:DP:SP:GQ  0/1:215,0,221:28:4:99   0/1:219,0,208:24:0:99
      #10      68575   .       C       T       999     .       DP=72;AF1=1;AC1=4;DP4=4,6,24,25;MQ=29;FQ=-60.5;PV4=0.73,1,1,1   GT:PL:DP:SP:GQ  1/1:235,29,0:31:4:63    1/1:239,41,0:28:2:75

      format = cols[8].split(':')
      
      if firstline == 1 
        i = 0
        format.each do |k|
          if k == "GT"
            fields[:gt] = i
          elsif k == "PL" 
            fields[:pl] = i
          elsif k == "SP"
            fields[:sp] = i
          end
          i+=1
        end
      end

      firstline = 2
      

      gt_tumor,gt_norm =  cols[9].split(':'), cols[10].split(':')
      
      if settings["--normal"] != 0 
        gt_tumor,gt_norm = gt_norm, gt_tumor
      end
      
      inclusionflag = 0     
      exclusionarray = []
      if gt_norm[fields[:gt]] == "0/0" and gt_tumor[fields[:gt]] == "0/1"
        inclusionflag = 1
      end
   
      if inclusionflag == 1
        npl = gt_norm[fields[:pl]].split(',')
        tpl = gt_tumor[fields[:pl]].split(',')

        ## normal sample: het pl > cutoff; tumor sample: ref sample > cutoff
        if npl[1].to_i < settings["--minPL"] 
          exclusionarray << "normal_low-PL"
        end
        
        if  tpl[0].to_i < settings["--minPL"] 
          exclusionarray << "tumor_low-PL"
        end
        
        info.each do |item|
          k,v=item.split('=')[0..1]
          if k == "DP4" 
            v1=v.split(',')
            if v1[2].to_f < settings["--dp4"] or v1[3].to_f < settings["--dp4"] 
              exclusionarray << "lowAlleleCountonEitherStrain"
            end
          elsif k == "PV4" 
            v1=v.split(',')
            if v1[0].to_f < settings["--pv41"] or v1[1].to_f < settings["--pv42"] or v1[2].to_f < settings["--pv43"] or v1[3].to_f < settings["--pv44"]
              exclusionarray << "PV4_significant"
            end
          elsif k == "CLR" 
            if v.to_f < settings["--clr"]
              exclusionarray << "lowCLR"
            end
          end
        end
      
      end

      
      if inclusionflag == 1 and exclusionarray.size == 0
        o.puts line
      else
        if exclusionarray.size > 0
         #  $stderr.puts "#{cols[1..4].join("\t")}\t#{exclusionarray.join(";")}"
          cols[6] = "#{pass};#{exclusionarray.join(";")}"
        end
        
        e.puts cols.join("\t")
      end
    end
  end
  o.close
  e.close
end

def countSamples(vcf)
  n = 0
  File.new(vcf, 'r').each do |line|
    n += 1
    if line.match("^#CHROM") 
      cols = line.chomp.split(/\s+/)
      return cols.size - 9 
    elsif n > 1000
      return 0
    end
  end
end

def getopt
  
  opts = GetoptLong.new(
                        ["--vcf", "-v", GetoptLong::REQUIRED_ARGUMENT],
                        ["--dp4", "-d", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--pv41", "-p", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--pv42", "-q", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--pv43", "-r", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--pv44", "-s", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--clr", "-c", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--normal", "-n", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--minPL", "-t", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--help", "-h", GetoptLong::NO_ARGUMENT]
                        )
  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  if optHash.key?("--help") or !optHash.key?("--vcf") 
    $stderr.puts "Usage: ruby __.rb --vcf VCF --normal [1/2] [--dp4[2] --pv41[0.0001] --pv42[0.0001] --pv43[0.0001] --pv44[0.0001] --clr[20] --minPL[40]]"
    $stderr.puts "     options: "
    exit
  end
  return optHash
  
end


main()


