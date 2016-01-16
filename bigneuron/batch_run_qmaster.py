import os
import time

c=0
os.system('ls *.qsub>./jobs.txt')
f= open("./jobs.txt")
j=0
for file_qsub in f:
        j = j+1
        c = c+1
        if c == 500:
             print j
             print "sleep for half and hour"
             time.sleep(60*30)
             c=0
        #print "qsub "+file_qsub
        os.system("qsub "+file_qsub)

print "done"
print str(j)+" jobs submitted"


