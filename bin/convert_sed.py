import os
import sys

if len(sys.argv) > 1:
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    output_file2 = sys.argv[3]
else:
    print("No arguments provided.")
    exit()

one = 1.0
ok = True
there = False
ier = 0
sfound = 0
ul = '  '

if not os.path.exists(input_file):
    print(f"File {input_file} not found")
    exit()

with open(input_file, "r") as lu_in:
    with open(output_file, "w") as lu_out:
        with open(output_file2, "w") as file12:
            file12.write("freq. ,flux ,err_flux ,MJD_start ,MJD_end ,flag,catalog,reference\n")
            line = lu_in.readline()
            sfound, stringin, str2, rra, rdec, rtype , rt = line.split()
            line = lu_in.readline()
            line = lu_in.readline()
            line = lu_in.readline()
            while ok:
                string = lu_in.readline()
                if string == "":
                    break
                if string[1:2] != "=":
                    freq, flux, err_up, err_lo, mjdstart, mjdend, flag, catalog, *rest = string.split()
                    flux_err = abs((float(err_up) - float(err_lo))) / 2.0
                    ll = len(rest) 
                    reference = ''
                    for i in range(0,ll):
                       reference = reference+' '+rest[i]
                    reference = reference.replace(",", " ")
                    if (float(flux) != 0.0) or (float(err_up) != 0.0):
                        if float(flux) < 0.0:
                            flux = -float(flux)
                        if flag == "UL":
                           ul = "UL"
                           if flux_err > 0:
                               flux = 2.0 * flux_err
                           flux_err = 0.0
                           lu_out.write(f"{float(rra):.5f} | {float(rdec):.5f} | {float(freq):.5e} | {one:.5e} | {float(flux):.5e} | {flux_err:.3e} | {float(mjdstart):.4f} | {float(mjdend):.4f} | {ul} | \n")
                           file12.write(f" {float(freq):.3e}, {float(flux):.3e}, {flux_err:.3e},  {float(mjdstart):.4f},  {float(mjdend):.4f},{ul},{catalog}     , {reference} \n")
                        else:
                           ul = " "
                           if (flux_err > float(flux)):
                               flux_err = min(float(err_up) - float(flux), float(flux) - float(err_lo))
                           lu_out.write(f"{float(rra):.5f} | {float(rdec):.5f} | {float(freq):.5e} | {one:.5e} | {float(flux):.5e} | {flux_err:.3e} | {float(mjdstart):.4f} | {float(mjdend):.4f} | {ul} | \n")
                           file12.write(f" {float(freq):.3e}, {float(flux):.3e}, {flux_err:.3e},  {float(mjdstart):.4f},  {float(mjdend):.4f}, {ul},{catalog}    , {reference} \n")
