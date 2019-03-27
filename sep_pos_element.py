import numpy as np
import os
import sys
import math
import shutil

#######################################
def get_atom_and_num(file_target):
    ele = []
    with open(file_target, "r") as p:
        find = False
        first_time = True
        for i, line in enumerate(p):
            if i == 5:
                ele += line.split()
    num = []
    with open(file_target, "r") as p:
        find = False
        first_time = True
        for i, line in enumerate(p):
            if i == 6:
                num += line.split()
    number_list = []
    for i in range(0,len(num)):
        number_list.append(int(num[i]))
    tot = sum(number_list)

    full_list = [None]*(3*tot)
    number_list.insert(0,0)
    number_list.append(0)
    for i in range(0, len(ele)):
        for x in range(3*sum(number_list[:(i+1)]),3*sum(number_list[:(i+2)])):
            full_list[x] = ele[i]

    return full_list
#######################################
def get_T_num(file_target,ele_list):
    counter = 0
    str1 = []
    ele_with_T = []
    with open(file_target, "r") as p:
        find = False
        first_time = True
        for i, line in enumerate(p):
            if ("Direct" in line) or ("Cartesian" in line):
                number = i
                find = True

            if find:
               if i > number:
                   b = line.split()[3:6]
                   if line.split()[3:6] != [None]:
                       str1 += b
    for i in range(0,len(str1)):
        if str1[i] == "T":
            counter += 1
    for i in range(0,len(str1)):
        if str1[i] == "T":
            ele_with_T.append(ele_list[i])


    return counter, str1, ele_with_T

######################################
def get_new_num_list(the_list):
    return_list = []
    num_list = []
    for i in range(1,len(the_list)):
        if the_list[i] != the_list[i-1]:
            return_list.append(i)
    return_list.insert(0,0)
    for i in range(0,len(return_list)-1):
        num_list.append(return_list[i+1]-return_list[i])
    num_list.append(len(the_list)-return_list[-1])
    return num_list
    
####################################
def T_num_in_folder(the_t, the_inp):

    folder = int(the_t/the_inp)
    the_list = [folder] * the_inp
    if the_t % the_inp != 0 :
        rest = the_t % the_inp
        for i in range(0, rest):
        
            the_list[i] += 1

    return the_list
###################################
def get_T_where(the_list):
    a = the_list
    return_list = []
    for i in range(0,len(the_list)):
        if the_list[i] == "T":
            return_list.append(i)
    return return_list
##################################

def T_sep_loc(sep,loc):
    return_list = [None] * len(new_num_list)
    ori = sep
    sep.insert(0,0)  
    for i in range(0,len(ori)-1):
        if i != len(ori)-2:
            return_list[i] = loc[sum(sep[0:i+1]):sum(sep[0:i+2])]
        if i ==len(ori)-2:
            return_list[i] = loc[sum(sep[0:i+1]):sum(sep)]
    sep.remove(0)
    return return_list
#################################
def get_coord(file_target):
    str1 = []
    with open(file_target, "r") as p:
        find = False
        first_time = True
        for i, line in enumerate(p):
            if ("Direct" in line) or ("Cartesian" in line):
                number = i
                find = True

            if find:
               if (number+int(len(TFlist)/3)+1>i)and(i > number):
                   str1 += line.split()[0:3]
    strr = np.array(str1).reshape((int(len(TFlist)/3)),3).tolist()

    return strr

######################################
def get_front(the_file):
    with open(the_file, "r") as p:
        find = True
        str3 = ""
        for i, line in enumerate(p):

            if find:
                str3 += line
                if ("Direct" in line) or ("Cartesian" in line):
                    find = False
    return str3

######################################
def new_pos(name,the_list): 
    F_list = ["F"] * len(TFlist)
    for a in the_list:
        F_list[int(a)] = "T"
    F_ = np.array(F_list).reshape((int(len(TFlist)/3)),3).tolist()
    f = open(name, 'w')
    f.write(front)
    str2 = ""
   
    for x in range(0, len(coord_list)):
        b = []
        b += coord_list[x]
        b += F_[x]
        str1 = ""
        for y in b:
            str1 += "\t" + y
        str2 += str1 + "\n"
    f.write(str2)
    f.close()


#################################
def change_INCAR(file,a,b,c,d,e,f,g,h):
    to_change = {a:b,c:d,e:f,g:h}
    incar_list = []
    with open(file, "r") as incar:
        for i, line in enumerate(incar):
            to_delete = False
            for x in to_change.keys():
                if x in line:
                    to_delete = True
            if not to_delete:
                incar_list.append(line)

    with open(file, "w") as incar_w:
        for x in incar_list:
            incar_w.write(x)
        for  y in to_change.keys():
            incar_w.write(to_change[y])
####################################################

def ori_folder(name):
    os.mkdir(name)
    shutil.copy2("INCAR",name)

    shutil.copy2("POSCAR",name)
    os.chdir(name)
    change_INCAR("INCAR","ISTART", "\nISTART = 0 \n","NSW","\nNSW = 1\n","LWAVE","\nLWAVE = .TRUE.\n","IBRION","\nIBRION = 1\n")

    os.symlink("../KPOINTS","KPOINTS")
    os.symlink("../POTCAR","POTCAR")
    os.chdir("../")
#################################
def all_new_pos(the_list):
    for i in range(0,len(the_list)):
        os.mkdir("0"+str(i))
        new_pos("00"+str(i),the_list[i])
        shutil.move("00"+str(i),"0"+str(i))
        shutil.copy2("INCAR","0"+str(i))
    for i in range(0,len(the_list)):
        os.chdir("0"+str(i))
        os.symlink("../KPOINTS","KPOINTS")
        os.symlink("../POTCAR","POTCAR")
        os.symlink("../original/WAVECAR","WAVECAR")
        change_INCAR("INCAR","ISTART", "\nISTART = 1 \n","NSW","\nNSW = 200\n","LWAVE","\nLWAVE = .FALSE.\n","IBRION","\nIBRION = 5\n")
        for y in os.listdir("."):
            if y == "00"+str(i):
                os.rename("00"+str(i),"POSCAR")
        os.chdir("../")
#################################

expend_ele_list = get_atom_and_num("POSCAR")
Tnum, TFlist, ele_with_T = get_T_num("POSCAR",expend_ele_list)
new_num_list = get_new_num_list(ele_with_T)
T_loc = get_T_where(TFlist)
sep_loc = T_sep_loc(new_num_list,T_loc)
coord_list = get_coord("POSCAR")
front = get_front("POSCAR")
ori_folder("original")
all_new_pos(sep_loc)
