# DESCRIPTION: python merge.py -a filename1 -b filename2 -c filename3 --file4=filename4

import getopt 
import sys 
 
config = { 
    "file1":"file1.txt", 
    "file2":"file2.txt", 
    "file3":"file3.txt", 
    "file4":"file4.txt", 
     
}
opts, args = getopt.getopt(sys.argv[1:], 'ha:b:c:d:',  
      [ 
        'file1=',  
        'file2=',  
        'file3=',
	'file4='
        ] 
      )

for option, value in opts: 
    if  option in ["-h","--help"]: 
        print """ 
        default help
        """ 
    elif option in ['-a','--file1']: 
        config["file1"] = value
    elif option in ['-b','--file2']: 
        config["file2"] = value
    elif option in ['-c','--file3']: 
        config["file3"] = value
    elif option in ['-d','--file4']: 
        config["file4"] = value

dict_a = {};
dict_b = {};
dict_c = {};
dict_key = {};
f = open(config["file1"])
line = f.readline()
while 1:
    line = line.strip()
    current_line = line.split('\t')
    key =current_line[0]+','+current_line[1]
    key_val =current_line[0]+current_line[1]
    dict_a[key] = [current_line[2],current_line[3]]
    dict_key[key_val] = key
    line = f.readline()
    if not line:
        break
    pass
f.close()

f = open(config["file2"])
line = f.readline()
while 1:
    line = line.strip()
    current_line = line.split('\t')
    key =current_line[0]+','+current_line[1]
    key_val =current_line[0]+current_line[1]
    dict_b[key] = [current_line[2],current_line[3]]

    dict_key[key_val] = key
    line = f.readline()
    if not line:
        break
    pass
f.close()

f = open(config["file3"])
line = f.readline()
while 1:
    line = line.strip()
    current_line = line.split('\t')
    key =current_line[0]+','+current_line[1]
    key_val =current_line[0]+current_line[1]
    dict_c[key] = [current_line[2],current_line[3]]

    dict_key[key_val] = key
    line = f.readline()
    if not line:
        break
    pass
f.close()
sort = sorted(dict_key.keys())
f = open(config["file4"], 'w')
for key_val in sort:
    key_s = dict_key[key_val].strip()
    key_array = key_s.split(',')
    if not dict_key[key_val] in dict_a:
        dict_a[dict_key[key_val]]=['0','0']
    if not dict_key[key_val] in dict_b:
        dict_b[dict_key[key_val]]=['0','0']
    if not dict_key[key_val] in dict_c:
        dict_c[dict_key[key_val]]=['0','0']        
#    print (key_array[0]+"\t"+key_array[1]+"\t"+dict_a[dict_key[key_val]][0]+"\t"+dict_a[dict_key[key_val]][1]+"\t"+dict_b[dict_key[key_val]][0]+"\t"+dict_b[dict_key[key_val]][1]+"\t"+dict_c[dict_key[key_val]][0]+"\t"+dict_c[dict_key[key_val]][1])
    f.write(key_array[0]+"\t"+key_array[1]+"\t"+dict_a[dict_key[key_val]][0]+"\t"+dict_a[dict_key[key_val]][1]+"\t"+dict_b[dict_key[key_val]][0]+"\t"+dict_b[dict_key[key_val]][1]+"\t"+dict_c[dict_key[key_val]][0]+"\t"+dict_c[dict_key[key_val]][1]+"\n")
f.close()
