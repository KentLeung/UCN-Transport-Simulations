import sys

num = int(sys.argv[2])

header1 = "connex: The first region should be region 0, and a special-handling code must be used to specify how to treat its cut-plane.\n"
header2 = "-----------------------------------------------------------------------------------------------------------------------------\n"
header3 = "Region#     Connects through its own cut-plane to:     AND also to:     AND also to:     AND also to:     AND also to:     AND also to:\n"
line1 = "0           0                                             1/\n"
output = [header1,header2,header3,line1]

for i in range(1,min(10,num)):
    output.append("{:d}{:12d}{:46d}/\n".format(i,i-1,i+1))

if num >= 10:
    for i in range(10,num):
        output.append("{:d}{:11d}{:46d}/\n".format(i,i-1,i+1))
    output.append("{:d}{:11d}/\n".format(num,num-1))
else:
    output.append("{:d}{:12./d}/\n".format(num,num-1))


file = open(sys.argv[1], "w")
file.writelines(output)
file.close()
