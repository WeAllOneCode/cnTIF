#create the interaction matrix from a cluster file
# and a summary contact as described below.
#example of cluster file in the input file: (gzip)
#The 7th column is interaction count
#chr1	3005312	3006423	chr1	4148852	4149854	1
#chr1	3004501	3005612	chr1	51631220	51632231	3
#chr1	3007277	3008278	chr1	5159978	5161089	19
#Output is a matrix, index & border files.
#filter out lower than mincount reads 

##########################################
import re
import sys
import math
import time
import gzip

#######################################
def outPutIndex(binning,mergeNum,outIndex):
    out=open("%s_i%d.ind"%(outIndex,mergeNum),'w')
    count=1
    out.write("INDEX%3d CHR      Start      End\n"%(mergeNum))
    for i in range(len(binning)):
        for j in range(len(binning[i])-1):
            out.write("%5d  %6s %10s %10s\n"%(count,get_chr_name(i+1),binning[i][j],binning[i][j+1]))
            count=count+1

    out.close()

    out=open("%s_i%d.border"%(outIndex,mergeNum),'w')
    count=1
    out.write("CHR  Start_Bin End_Bin Length\n")
    for i in range(len(binning)):
        out.write("%5s %6s %6s %6s\n"%(get_chr_name(i+1),count,count+len(binning[i])-2,len(binning[i])-1))
        count=count+len(binning[i])-1

    out.close()

def assign_matrix(matrix,chr_bead,bead_len):
    chrr=1
    count=1
    for m in range(1,bead_len+1):
        matrix[0][m]=matrix[m][0]=count
        count=count+1
        if m==chr_bead[chrr]:
            chrr=chrr+1
            count=1

def get_chr_name(num):
    if num==23:
        return "chrX"
    elif num==24:
        return "chrY"
    elif num<23:
        return "chr"+str(num)
    else:
        print "Wrong Chromosome Name"
        sys.exit()

def get_chr_num(ch):
    try:
        chnum = int(ch)
    except:
        if ch == 'X':
            chnum = 23
        elif ch == 'Y':
            chnum = 24
        else:
            chnum = 25
    return chnum

def get_bead_id(ch,p,binning,chr_bead):
    #(ch-1) is used because the python list begins with 0
    binning_chr=binning[ch-1] #binned sites
    pos = int(p)
    for i in range(1,len(binning_chr)):
        if pos<=binning_chr[i] and pos>binning_chr[i-1]:
            #print 'Got bin:',i, binning_chr[i-1],pos,binning_chr[i],ch
            binning_pos=i
            break
    else:
        print "Can not assign bead ID! Please check the following:"
        print 'Last index & bin evaluated for ch & position:',i,binning_chr[-1], ch, pos
        sys.exit()
    return binning_pos+chr_bead[ch-1]

#######################################
file=sys.argv[1]  #gz file
runID=sys.argv[2] # just a label
bpbin=int(sys.argv[3]) # bin size
mincount = 0 # minimum iPET

short_kb = 8000 #20kb cutoff filter
exclude_chr = ["Y","M"] #not will be in matrix; input summary file does not contain 'chr'-string
#######################################

chr_length=[249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560]

myBinning = [] #all borders of chromosomes


for cl in chr_length:
    useq = [ j*bpbin for j in range(2+cl/bpbin)]
    myBinning.append(useq)

outPutIndex(myBinning,bpbin,runID) #just print out the ind and border files

print "Wrote indexing %s_i%d"%(runID, bpbin)
#sys.exit()


chr_bead=[0]
for i in range(len(myBinning)):
    chr_bead.append(chr_bead[i]+len(myBinning[i])-1)
bead_len=chr_bead[-1]

t1=time.time()

matrix=[[0]*(bead_len+1) for i in range(bead_len+1)]
assign_matrix(matrix,chr_bead,bead_len)
#print matrix[0] #this is the header of matrix


print 'Reading %s'%(file)

if file[-2:]=='gz':
    f=gzip.open(file,'rb')
else:
    f=open(file,'r')

line=f.readline()

counts=0
icounts=0
total_without_20kb=0
intra_without_20kb=0
interchr=0
lineNum=0
mitoch=0

includeChrs = [str(i) for i in range(1,23)] #adapting for itx file that contains unused chrs
includeChrs.append('X')

while line:
    body=filter(None,re.split(' |\t|\n',line))
    ch1 = re.sub("chr","",body[0])
    ch2 = re.sub("chr","",body[3])
    itx=int(body[6])
    if (itx > mincount and ch1 in includeChrs and ch2 in includeChrs):
        mid1 = 0.5*(int(body[1]) + int(body[2]) )
        mid2 = 0.5*(int(body[4]) + int(body[5]) )
        if ch1==ch2 and abs(mid2-mid1)<=short_kb:
            pass
        else:
            if ch1==ch2:
                intra_without_20kb=intra_without_20kb+itx
            else:
                interchr=interchr+itx
            total_without_20kb=total_without_20kb+itx
            try:
                bead1=get_bead_id(get_chr_num(ch1),mid1,myBinning,chr_bead)
                bead2=get_bead_id(get_chr_num(ch2),mid2,myBinning,chr_bead)
            except:
                print 'WRONG input !!!!!!'
                print line
                sys.exit()
            matrix[bead1][bead2]=matrix[bead1][bead2]+itx
    else:
        mitoch=mitoch+itx
    counts=counts+1 
    icounts=icounts+itx 
    line=f.readline()
    lineNum=lineNum+1
    if lineNum%1000000==0:
        t3=time.time()
        print "%d M in %.1f s"%(lineNum/1000000,t3-t1)
f.close()
print "Total PET lines & interaction counts: %d %d"%(counts, icounts)
print "Total itx passed filter:",total_without_20kb
print "Intra without self-ligation:",intra_without_20kb
print "Inter ",interchr

for i in range(1,len(matrix)):
    for j in range(i+1,len(matrix[i])):
        matrix[i][j]=matrix[j][i]+matrix[i][j]
        matrix[j][i]=matrix[i][j]

summ=0
for i in range(1,len(matrix)):
    summ=summ+sum(matrix[i][1:])+matrix[i][i] #this is double counting way
print "checking total itxs in matrix= %d; #of other chrs:%d"%(summ/2,mitoch)

out=open('%s_i%d_rawPET%d.mat'%(runID, bpbin, mincount),'w')
for j in range(len(matrix)):
    for k in range(len(matrix[j])):
        if j*k==0:
            out.write('%9s'%(matrix[j][k]))
        else:
            if matrix[j][k]==0:
                out.write('%9s'%(matrix[j][k]))
            else:
                out.write('%9s'%(matrix[j][k]))
    out.write('\n')
out.close()

t3=time.time()
print 'Total time spent: %.2f sec'%(t3-t1)
sys.exit()
#--------------------------------------------------------------------

