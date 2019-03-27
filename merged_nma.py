import numpy
import math
import sys
import os
import shutil
##############
# constats
h_IN_Js                = 6.626070040e-34
eV_IN_J                = 1.60218e-19
angstom_IN_m           = 1.0e-10 
amu_IN_kg              = 1.66054e-27
speed_of_light_m_per_s = 299792458.0

expression_fix= " ".join(",".join(sys.argv[1:]).split(",")).split()

##############
def get_all_out(directory):

    return_list = []
    the_dir = []
    for x in os.listdir(directory):
        if x[:1] == "0":
            the_dir.append(x)

    for y in range(0,len(the_dir)):
        for z in os.listdir(directory+"/"+the_dir[y]):
            if z == "OUTCAR":
                os.chdir(directory+"/"+the_dir[y])
                return_list.append(os.getcwd()+"/"+z)
                os.chdir("../")
    return return_list

##############
def is_nth_coordinate(atom_number, xyz):
    return_list = []
    for i in range (0,len(xyz)):
        if xyz[i] == "x":
            return_list.append(( (atom_number-1)*3)   )
        elif xyz[i] == "y":
            return_list.append(( (atom_number-1)*3) +1)
        elif xyz[i] == "z":
            return_list.append(( (atom_number-1)*3) +2)

    return return_list
##############
def string_is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False
##############
def interpret_from_to_statement(s):
    start = s.split("-")[0]
    end   = s.split("-")[1]
    return_coords = []

    s_n, s_c = separate_atom_number_and_xyz(start)   # Start_atom_Number , Start_Coordinate
    e_n, e_c = separate_atom_number_and_xyz(end)   # end_ , end_

    if s_c != e_c:
        print ("\n\nCannot understand Expression {0}.".format(s))
        print ("""Examples for allowed in put type:
1,2,3, 4-23  ---> fix xyz of the atoms
3x,4yx ---> fix x coordinate of 3, x and y of 4
3x-6x  ---> fix x coordiante of all atoms 3 to 6""")

        quit()
    elif (s_n < 0) or (e_n < 0):
        print ("\n\nCannot understand Expression {0}.".format(s))
        print ("""Examples for allowed in put type: 
1,2,3, 4-23  ---> fix xyz of the atoms
3x,4yx ---> fix x coordinate of 3, x and y of 4
3x-6x  ---> fix x coordiante of all atoms 3 to 6""")
        quit()
    elif s_n < e_n:
        for x in range(s_n, e_n+1):
            return_coords += is_nth_coordinate(x, s_c)
    else:
        for x in range(e_n, s_n+1):
            return_coords += is_nth_coordinate(x, s_c)
    
    return return_coords

##############
def separate_atom_number_and_xyz(s):
    for x in range(len(s), 0 , -1):
        if string_is_int(s[0:x]):
            if x== len(s):
                return int(s[0:x]), "xyz"
            else:
                return int(s[0:x]), s[x:]
        else:
            if not s[x-1] in ["x","y","z"]:
                print ("\n\nCannot understand Expression {0}.".format(s))
                print ("""Examples for allowed in put type: 
1,2,3, 4-23  ---> fix xyz of the atoms
3x,4yx ---> fix x coordinate of 3, x and y of 4
3x-6x  ---> fix x coordiante of all atoms 3 to 6""")
                quit()
##############

def get_fixed_coordinates(expression_fix):
    coords_to_fix = []
    for x in expression_fix:
        if "-" in x:
            coords_to_fix += interpret_from_to_statement(x)
        else:
            number, xyz = separate_atom_number_and_xyz(x)
            coords_to_fix += is_nth_coordinate(number, xyz)

    return coords_to_fix

##############


def str_list2int_list(str_list):
    return_list = []
    for x in str_list:
        return_list.append(int(x))

    return return_list
##############
def str_list2float_list(str_list):
    return_list = []
    for x in str_list:
        return_list.append(float(x))

    return return_list
##############
def number_list2str_list(number_list, digits):
    return_list = []
    form = "." + str(digits) + "f"
    for x in number_list:
        return_list.append(  format(x, form ))

    return return_list
##############
def get_mass_and_force(out_list):
    ion_found = False
    mass_found = False
    mass_au_list = []
    start_read_force = False
    current_force = []
    all_forces = []
    for num in range(0,len(out_list)):
        mass_au = []
        numbers = []
        with open(out_list[num],"r") as o:
            for i,line in enumerate(o):
                if not (ion_found and mass_found):
                    
                    if "ions per type" in line:
                        numbers+=str_list2int_list(line.split("=")[1].split())
                        ion_found = True
                    elif ("POMASS" in line) and ("ZVAL" in line):
                        mass_au.append(float((line.split("=")[1]).replace("; ZVAL","")))    
                    if len(mass_au) == len(numbers):
                        
                       # if ion_found and mass_found:
                        if num == 0:
                            for j in range(0,len(numbers)):
                                for k in range(0, numbers[j]*3):
                                    mass_au_list.append(mass_au[j])
                        mass_found = True
                else:
                    if start_read_force:
                        if not ("----------------------------" in line) and not ("drift"  in line):
                            current_line = line.split()
                            current_force += [[current_line[0], current_line[3] ],]
                            current_force += [[current_line[1], current_line[4] ],]
                            current_force += [[current_line[2], current_line[5] ],]

                    if "TOTAL-FORCE" in line:
                        start_read_force = True
                        current_force= []
                    if "total drift" in line:
                        start_read_force = False
                        all_forces.append(current_force)

    return mass_au_list, all_forces   # list of  atomic  mass in au, all masses appear 3 times as, one each direction ( xyz)
                                      # all forces is a list of list of list, mittle list are the forces calculated for displacements of differnet coordinates,
                                      # element 1 of the inner list is the position of the n-th atomic cartesian coordiante, and the second element is 
                                      # force along this direction
################
###################################
def analyze_force_list(force_list):   #here the programm first trieds to identify the equilibrium geometry, when reading from multiple OUTCAR files, the equilibrium geometric 
                                      #not necessary the first entry
    changed_coords_unsort = []
    fixed_coord = []

    for x in force_list:
        current = []
        for j in range(0, len(force_list[0])):
            if force_list[0][j][0] != x[j][0]:
                current.append(j)                # this is a list that says: x line in force_list is a displacement of the j-th coordinate, if equilibrium geometry is not the first i
                                                 # there may be two entries in current 
        changed_coords_unsort.append(current)    # list of list with current 



#    print changed_coords_unsort

    if eq_is_first(changed_coords_unsort):
        for count in range(0, len(force_list[0])):
            if not [count, ] in changed_coords_unsort:
                fixed_coord.append(count)
        return changed_coords_unsort, fixed_coord
    else:
        change_one = []
        change_zero = []
        minima = []
        changed_coords_sort = []
        for i in range(0,len(changed_coords_unsort)):
            if len(changed_coords_unsort[i]) == 1 :
                change_one.append(i)
            elif len(changed_coords_unsort[i]) == 0 :
                change_zero.append(i)
            elif len(changed_coords_unsort[i]) == 2:
                has_two_different = i
        to_remove = changed_coords_unsort[change_one[0]][0]
        for k in change_one:
            current = []
            for  l in range(0, len(force_list[k])):
                if force_list[has_two_different][l][0] != force_list[k][l][0]:
                    current.append(l)
            if len(current) == 1:
                minima.append(k)

        for m in range(0,len(changed_coords_unsort)):
            if m in change_one:
                if m in minima:
                    changed_coords_unsort[m] = []
            elif m in change_zero:
                changed_coords_unsort[m] = [to_remove]
            else:
                new = remove_from_list(changed_coords_unsort[m], to_remove)
                changed_coords_unsort[m] = new
    for count in range(0, len(force_list[0])):
        if not [count, ] in changed_coords_unsort:
            fixed_coord.append(count)

    return changed_coords_unsort, fixed_coord    #changed_coords_unsort is list of list, it shows which coordinate has been displaced at which geommetry. 
                                                 #the inner lists should be list with only one element 
                                                 #fixed_coord is a list with all fixed coordinates
###############
def eq_is_first(changed_coord):
    for x in changed_coord:
        if len(x) > 1:
            return False
    return True

###############
###############
def remove_from_list(l, to_remove):
    return_list  = []
    for x in l:
        if x != to_remove:
            return_list.append(x)
    return return_list
##############
def which_pertubation_force(l, pertubation_list):
    return_list = []
    for x in range(0,len(l)):
        if l[x] == pertubation_list:
            return_list.append(x)
    return return_list
###############
def get_ave_in_list(x_list,y_list):
    same_list = []
    for i in range(0, len(x_list)):
        rem_list = []
        first_part = x_list[:i]
        last_part = x_list[(i+1):]
        rem_list += first_part
        rem_list += last_part
        for z in range(0, len(rem_list)):
            if x_list[i] == rem_list[z]:
                same_list.append(x_list[i])

    a = same_list
    seen = {}
    dupes = []
    for x in a:
        if x not in seen:
            seen[x] = 1
        else:
            if seen[x] == 1:
                dupes.append(x)
            seen[x] += 1

    new_y_list = []
    new_x_list = []
    index_list = []
    for i in range(0,len(dupes)):
        counter = 0
        tot_x = 0
        tot_y = 0
        for x in range(0,len(x_list)):
            if x_list[x] == dupes[i]:
                tot_x += float(x_list[x])
                tot_y += float(y_list[x])
                counter += 1
                index_list.append(x)
        ave_x = tot_x / counter
        ave_y = tot_y / counter

        new_y_list.append(str(ave_y))
        new_x_list.append(str(ave_x))

    for a in range(0, len(y_list)):
        if a not in index_list:
            new_y_list.append(y_list[a])
            new_x_list.append(x_list[a])    


    return new_x_list,new_y_list
###############
def get_second_deriv(changed_coord, forces_list, fixed, vari):

    return_deriv = []

    x_val = []
    y_val = []

    place_perturb_data = which_pertubation_force(changed_coord,[]) + which_pertubation_force(changed_coord,[vari])

#    print place_perturb_data     
    perturb_data = []
    for x in place_perturb_data:
        perturb_data.append(forces_list[x])



    for y in range(0, len(perturb_data[0])):
        if not  y in fixed:
            x_val = []
            y_val = []
            for z in range(0, len(perturb_data)):
                x_val.append(perturb_data[z][vari][0])
                y_val.append(perturb_data[z][y][1])
            x_val_new, y_val_new = get_ave_in_list(x_val,y_val)
            x_arr = numpy.array(str_list2float_list(x_val_new))
            y_arr = numpy.array(str_list2float_list(y_val_new))
            grad, y_intercept = numpy.polyfit(x_arr,y_arr,1)
            return_deriv.append(grad)

    return return_deriv
###############

def get_hessian(changed_coord, force_list, fixed_coord):
    hessian = []
    for x in range(0, len(force_list[0])):
        if not x in  fixed_coord:
            a = get_second_deriv(changed_coord, force_list, fixed_coord, x)
            hessian.append(a)


    make_matrix_sym(hessian)

    return numpy.matrix(hessian)
###############
def make_matrix_sym(M):
# if M_ij and M_ji are not the same, replace both values with the average value
    for i in range(0,len(M)):
        for j in range(0,len(M[0])):
            if  i!=j:
                if M[i][j] !=  M[j][i]:
                    average =  (M[i][j] + M[j][i])/2
                    M[i][j] = average
                    M[j][i] = average

###############
def one_over_sqrt_mass_list_no_fixed(mass_all, fixed):
    return_list = []
    for i in range(0, len(mass_all)):
        if not i in fixed:
            return_list.append(1.0/ math.sqrt(mass_all[i]))

    return numpy.matrix(return_list)
###############
def mass_weight(hessian, one_over_sqrt_mass):
    return numpy.diag(one_over_sqrt_mass.A1) * ( hessian * numpy.diag(one_over_sqrt_mass.A1))
###############
def get_order(float_list):
    sort_l = []
    for x in range(0, len(float_list)):
        sort_l.append([x,float_list[x]])

    return_list = []
    for y in  sorted(sort_l, key=lambda entry: entry[1]):
        return_list.append(y[0])

    return return_list
###############
def get_freq_and_eigenvec(hessian, mass_au, fixed):

    perturbed_1_over_sqrt_mass_au = one_over_sqrt_mass_list_no_fixed(mass_au, fixed)

    mass_weighted_hessian = mass_weight(hessian, perturbed_1_over_sqrt_mass_au)

#    for x in hessian:
#        print x.tolist()
#
#    print 
#    for x in mass_weighted_hessian:
#        print x.tolist() 

    eigen_val, mass_weight_eigen_vec = numpy.linalg.eig(mass_weighted_hessian)

    eigen_val_list = eigen_val.tolist()
    order_of_eigen_val =   get_order(eigen_val_list)

#    cart_eigen_vec = numpy.diag(perturbed_1_over_sqrt_mass_au.A1) * mass_weight_eigen_vec
    cart_eigen_vec = mass_weight_eigen_vec 
    cart_eigen_vec_no_mass = numpy.diag(perturbed_1_over_sqrt_mass_au.A1) * mass_weight_eigen_vec 
 
    extended_modes = add_fixed_coord(cart_eigen_vec.T.tolist(), fixed)
    extended_modes_no_mass = add_fixed_coord(cart_eigen_vec_no_mass.T.tolist(), fixed)

    sorted_eigenval = []
    sorted_mode =     []
    sorted_mode_no_mass = []
    for i in order_of_eigen_val:
        sorted_eigenval.append(eigen_val_list[i]) 
        sorted_mode.append(extended_modes[i])  
        sorted_mode_no_mass.append(extended_modes_no_mass[i])
 
    return sorted_eigenval, sorted_mode, sorted_mode_no_mass

###############
def add_fixed_coord(cart_eigen_vec_list, fixed):
    extended_modes = []
    for x in cart_eigen_vec_list:
        counter_mode = 0
        current_extended_mode = []
        for y in range(0, len(x)+len(fixed)):
            if y in fixed:
                current_extended_mode.append(0.0)
            else:
                current_extended_mode.append(x[counter_mode])
                counter_mode += 1
        atom_wise_mode = []
        count = 0
        atom = " "
        #print current_extended_mode
        for xx in current_extended_mode:
            blank = " "*(5-len(str(xx).split(".")[0]))
            atom +=  blank + "{0:.6f}".format(xx)
            count += 1
            if count == 3 :
                atom_wise_mode.append(atom)
                count = 0
                atom = " "

        extended_modes.append(atom_wise_mode)

    return extended_modes

###############

def get_coord(changed_coord, force_list):

    return_list = []

    for x in range(0, len(changed_coord)):
        if changed_coord[x] == []:
            for i in force_list[x]:
                return_list.append(i[0])

            return return_list   #this return all cartesian coordinates of the equilibrium geometry 

###############
def three_entry_in_line(str_list):
    return_list = []
    counter = 0
    current_line = "    "
    for i in range(0, len(str_list)):
        space = " " * (3-len(str_list[i].split(".")[0]))
        zeros = "0" * (6-len(str_list[i].split(".")[1]))

        current_line += space + str_list[i] +zeros  
        counter += 1 
        if counter == 3 :
            return_list.append(current_line)
            current_line = "    "
            counter = 0

    return return_list 
###############
def three_entry_in_line_uniform_string(str_list):
    return_list = []
    counter = 0
    current_line = "    "
    for i in range(0, len(str_list)):
        space = " " * (3-len(str_list[i]))

        current_line += space + str_list[i]
        counter += 1
        if counter == 3 :
            return_list.append(current_line)
            current_line = "    "
            counter = 0

    return return_list
################
def make_output(output_file, eigen_val, atom_pos, eigen_mode, head):
    with open(output_file, "a") as out:
        out.write(head)
        out.write("\n Eigenvectors and eigenvalues of the dynamical matrix\n")
        out.write(" ----------------------------------------------------\n\n")
        for x in range(0, len(eigen_val)):
            blank  = "\n" + " " * (4-(len(str(1+x))))
            
            out.write(blank + str(1+x) + " " + freq_line(eigen_val[x]))
            out.write("       x         y         z           dx          dy          dz\n")

            #print atom_pos
            #print len(atom_pos)
            for y in range(0,len(atom_pos)):
                out.write(atom_pos[y] + eigen_mode[x][y] +"\n")

###############
def freq_line(freq_in_eV_overA2_amu):

    if freq_in_eV_overA2_amu < 0.0:
        head = "f  ="
        freq_in_eV_overA2_amu *= -1 
    else:
        head = "f/i="
        
    two_pi_THz_f = math.sqrt(freq_in_eV_overA2_amu*eV_IN_J/angstom_IN_m/angstom_IN_m/amu_IN_kg) * (1.0e-12)     # angular_freq
    two_pi_THz   = " {0:.6f}".format(two_pi_THz_f) + " 2PiTHz"
    blank = " " * (6-(len(two_pi_THz.split(".")[0])))
    two_pi_THz = blank + two_pi_THz


    THz_f =  two_pi_THz_f / (math.pi * 2.0  )                           # freq
    THz   = " {0:.6f}".format(THz_f) + " THz" 
    blank = " " * (5-(len(THz.split(".")[0])))
    THz = blank + THz

    meV_f = THz_f * 1.0e12 * h_IN_Js * 1000.0 /eV_IN_J
    meV   = " {0:.6f}".format(meV_f) + " meV"
    blank = " " * (6-(len(meV.split(".")[0])))
    meV = blank + meV

    wavenumber_f = THz_f *1.0e10 / speed_of_light_m_per_s 
    wavenumber = " {0:.6f}".format(wavenumber_f) + " cm-1"
    blank = " " * (5-(len(wavenumber.split(".")[0])))
    wavenumber = blank + wavenumber

#    print wavenumber_f/  math.sqrt(math.fabs(freq_in_eV_overA2_amu))
    return head + THz + two_pi_THz +  wavenumber + meV + "\n"

###############
def get_output_head(read_out, atom_pos):
    needed_lines = []
    lattice_read = False
    dir_and_reci_read = False
    vector = []
    with open(read_out, "r") as f:
        for i,line in enumerate(f):
            if i == 0 :
                needed_lines.append(line)
                needed_lines.append("INCAR:\n")
            if "TITEL" in line:
                needed_lines.append(line)
            if ("Lattice vectors:" in line) and not lattice_read:
                needed_lines.append(line)
                lattice_read = True
                vec_lines = [i+1, i+2, i+3, i+4]
            if lattice_read:
                if i in vec_lines:
                     needed_lines.append(line)
            if "ions per type" in line:
                needed_lines.append(line + "\n")
             

            if ("direct lattice vectors" in line) and not dir_and_reci_read:
                needed_lines.append("\n"+line)
                dir_and_reci_read = True
                dir_reci_lines = [i+1, i+2, i+3 ]
            if dir_and_reci_read:
                if i in dir_reci_lines:
                     needed_lines.append(line)
                     vector.append(  str_list2float_list(line.split()[0:3]))



    needed_lines.append("\n position of ions in fractional coordinates (direct lattice)\n")

                    
    cell_vec = numpy.matrix(vector)
    basis_trans = cell_vec.T.I 

    for atom in atom_pos:
        needed_lines.append("   " + " ".join(number_list2str_list(numpy.dot(basis_trans, numpy.array(str_list2float_list(atom.split()[0:3]))).tolist()[0], 8))+ "\n")
        
 
    return "".join(needed_lines)
  

###############
def element_list_form_output_head(head):
    elements=[]
    for x in head.split("\n"):
        if "TITEL" in x:
            elements.append(x.split()[3])
        if "ions per type" in x:
            occurence = x.split()[4:]

    elements_list = []
    for i in range(0, len(occurence)):
        for j in range(0, int(occurence[i])):
            elements_list.append(elements[i])

    return elements_list 


###############
def fix_atoms(list_to_fix, changed_coord, force_list, fixed_coord):
   
    list_of_list_to_fix = []
    for i in list_to_fix:
        list_of_list_to_fix.append([i,])

    new_changed_coord = []
    new_force_list    = []
 
    for x in range(0, len(changed_coord)):
        if changed_coord[x] in list_of_list_to_fix:
            pass
        else:
            new_changed_coord.append(changed_coord[x])
            new_force_list.append(force_list[x])


    relaxed_or_not = []
    for y in range(0, len(force_list[0])):
        if y in fixed_coord:
            relaxed_or_not.append("XX")
        elif y in list_to_fix:
            relaxed_or_not.append("# ")
        else:
            relaxed_or_not.append(". ")

    relaxed_or_not_lines = three_entry_in_line_uniform_string(relaxed_or_not)
 
        
    return new_changed_coord, new_force_list, sorted(list(set(list_to_fix + fixed_coord))), relaxed_or_not_lines
    

###############
OUT_list = get_all_out(os.getcwd())

mass_au_list, force_list =get_mass_and_force(OUT_list)


#temp_force_list = forces_list[1:]+[forces_list[0],]
#force_list = temp_force_list
changed_coord, fixed_coord = analyze_force_list(force_list)

atom_pos_output = three_entry_in_line(get_coord(changed_coord, force_list))   #  this is used for the output file, gets all coordinates 


list_to_fix = get_fixed_coordinates(expression_fix) 
changed_coord_new, force_list_new, fixed_coord_new, fixed_lines_print_out = fix_atoms(list_to_fix, changed_coord, force_list, fixed_coord)

#for x in range(0, len(force_list_new)):
#    print "\n\n"
#    for y in range(0,len(force_list_new[x])):
#        print  str(x+1) + "\t" + str(y+1) + "\t" + str(force_list_new[x][y][0]) + "\t" + str(force_list_new[x][y][1]) 

hessian = get_hessian(changed_coord_new, force_list_new, fixed_coord_new)


eigenval, eigenmode, eigenmode_nomass = get_freq_and_eigenvec(hessian, mass_au_list, fixed_coord_new)

#print eigenval

output_head = get_output_head(OUT_list[0], atom_pos_output)


element = element_list_form_output_head(output_head)

relax_lines  = "\n\n\n\n#    # = fixed, XX = not relaxed in original OUTCAR, . = relaxed\n\n\n"
for x in range(0, len(fixed_lines_print_out)):
    relax_lines += "## "+ str(x+1)+ "\t" +element[x]+"\t" +atom_pos_output[x] + fixed_lines_print_out[x]+ "\n"

relax_lines += "\n\n\n\n"
with open ("OUTCAR_norm","a") as o:
    o.write("--------------------------------------------------------------------------------------------------------------- "+ "\n")
    o.write("This program is written by Cheng-Chau Chiu and Yu-Chi Lee \n")
    o.write("If you use this program, please cite  Cheng-Chau Chiu, Yu-Chi Lee, and https://github.com/cccde/vasp_vibration_para\n")
    o.write("--------------------------------------------------------------------------------------------------------------- \n\n\n\n\n\n") 
make_output("OUTCAR_norm", eigenval, atom_pos_output, eigenmode, relax_lines + output_head)
#make_output("OUTCAR_NM", eigenval, atom_pos_output, eigenmode_nomass, relax_lines + output_head)

#element = element_list_form_output_head(output_head)

#relax_lines  = ""
#for x in range(0, len(fixed_lines_print_out)):
#    relax_lines += "## "+ str(x+1)+ "\t" +element[x]+"\t" +atom_pos_output[x] + fixed_lines_print_out[x]+ "\n"
