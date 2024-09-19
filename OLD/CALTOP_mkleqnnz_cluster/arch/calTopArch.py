from pyxdsm.XDSM import XDSM, OPT, SOLVER, FUNC, LEFT, GROUP

# Change `use_sfmath` to False to use computer modern
x = XDSM(use_sfmath=True)

x.add_system("main", OPT, r"\text{main.c}")
#x.add_system("openap", SOLVER, r"\text{OpenAP}")

# Define tier-one functions here
x.add_system("readinput", GROUP, r"\text{readinput.c}")
x.add_system("allocation", GROUP, r"\text{allocation.f}")
x.add_system("calinput", GROUP, r"\text{calinput.f}")
x.add_system("rho", GROUP, r"\text{rho.c}")
x.add_system("writeheading", GROUP, r"\text{writeheading.c}")
x.add_system("tiedcontact", GROUP, r"\text{tiedcontact.c}")
x.add_system("cascade", GROUP, r"\text{cascade.c}")
x.add_system("mastruct", GROUP, r"\text{mastruct.c}")
x.add_system("densityfilter", GROUP, r"\text{densityfilter.c}")
x.add_system("linstatic", GROUP, r"\text{linstatic.c}")
x.add_system("nonlingeo", GROUP, r"\text{nonlingeo.c}")
x.add_system("sensitivity", GROUP, r"\text{sensitivity.c}")
x.add_system("filterVector", GROUP, r"\text{filterVector.c}")

# Define second-tier functions here
#x.add_system("includefilename", FUNC, r"\text{includefilename.f}")
#x.add_system("getnewline", FUNC, r"\text{getnewline.f}")
#x.add_system("inputerror", FUNC, r"\text{inputerror.f}")
#x.add_system("nidentk", FUNC, r"\text{nidentk.f}")
#x.add_system("restartshort", FUNC, r"\text{restartshort.f}")
#x.add_system("userelements", FUNC, r"\text{userelements.f}")
#x.add_system("keystart", FUNC, r"\text{keystart.f}")
#x.add_system("splitline", FUNC, r"\text{splitline.f}")
#x.add_system("boundarys", FUNC, r"\text{boundarys.f}")
#x.add_system("buckles", FUNC, r"\text{buckles.f}")
#x.add_system("cloads", FUNC, r"\text{cloads.f}")
#x.add_system("dampings", FUNC, r"\text{dampings.f}")
#x.add_system("densitys", FUNC, r"\text{densitys.f}")
#x.add_system("depvars", FUNC, r"\text{depvars.f}")
#x.add_system("designvariables", FUNC, r"\text{designvariables.f}")
#x.add_system("dloads", FUNC, r"\text{dloads.f}")
#x.add_system("elastics", FUNC, r"\text{elastics.f}")
#x.add_system("noelfiles", FUNC, r"\text{noelfiles.f}")
#x.add_system("elprints", FUNC, r"\text{elprints.f}")
#x.add_system("noelsets", FUNC, r"\text{noelsets.f}")
#x.add_system("sectionprints", FUNC, r"\text{sectionprints.f}")
#x.add_system("materials", FUNC, r"\text{materials.f}")
#x.add_system("nodes", FUNC, r"\text{nodes.f}")
#x.add_system("nodeprints", FUNC, r"\text{nodeprints.f}")
#x.add_system("objective", FUNC, r"\text{objective.f}")
#x.add_system("outputs", FUNC, r"\text{outputs.f}")
#x.add_system("solidsections", FUNC, r"\text{solidsections.f}")
#x.add_system("statics", FUNC, r"\text{statics.f}")
#x.add_system("steps", FUNC, r"\text{steps.f}")
#x.add_system("surfaces", FUNC, r"\text{surfaces.f}")
#x.add_system("gen3delem", FUNC, r"\text{gen3delem.f}")
x.add_system("mafillsmmainVectorfilter", FUNC, r"\text{mafillsmmainVectorfilter.c}")

x.add_system("mafillsmmainVectorfiltermt", FUNC, r"\text{mafillsmmainVectorfiltermt.c}")
x.add_system("mafillsmvectorfilter", FUNC, r"\text{mafillsmvectorfilter.f}")
x.add_system("elementcpuload", FUNC, r"\text{elementcpuload.c}")

x.add_system("mafillsmmainse", FUNC, r"\text{mafillsmmainse.f}")
x.add_system("mafillsmse", FUNC, r"\text{mafillsmse.f}")
x.add_system("ec3dse", FUNC, r"\text{ec3dse.f}")

x.add_system("mafillsmmain", FUNC, r"\text{mafillsmmain.c}")
x.add_system("mafillsmas", FUNC, r"\text{mafillsmas.f}")
x.add_system("spooles", FUNC, r"\text{SPOOLES.c}")


# readinput.f -------------------------------#
#x.connect("readinput", "includefilename", r"\text{buff,..,includefn}")
#x.connect("readinput", "keystart", r"\text{ifreeinp,..,ikey}")
#x.connect("readinput", "splitline", r"\text{buff,..,n}")
#x.connect("readinput", "restartshort", r"\text{nset,..,nef}")
#---------------------------------------------#

# Allocation.f -------------------------------#
#x.connect("allocation", "getnewline", r"\text{inpc,..,ipoinpc}")
#x.connect("allocation", "inputerror", r"\text{inpc,..,iline}")
#x.connect("allocation", "nidentk", r"\text{iuel,..,four}")
#x.connect("allocation", "restartshort", r"\text{iuel,..,four}")
#x.connect("allocation", "userelements", r"\text{textpart,..,ipol}")
#---------------------------------------------#


# calinput.f -------------------------------#
#x.connect("calinput", "getnewline", r"\text{inpc,..,ipoinpc}")
#x.connect("calinput", "boundarys", r"\text{inpc,..,ier}")
#x.connect("calinput", "buckles", r"\text{inpc,..,ier}")
#x.connect("calinput", "cloads", r"\text{inpc,..,flag}")
#x.connect("calinput", "dampings", r"\text{inpc,..,nmat}")
#x.connect("calinput", "densitys", r"\text{inpc,..,ier}")
#x.connect("calinput", "depvars", r"\text{inpc,..,ier}")
#x.connect("calinput", "designvariables", r"\text{inpc,..,ier}")
#x.connect("calinput", "dloads", r"\text{inpc,..,ier}")
#x.connect("calinput", "elastics", r"\text{inpc,..,ier}")
#x.connect("calinput", "noelfiles", r"\text{inpc,..,ier}")
#x.connect("calinput", "elprints", r"\text{inpc,..,flag}")
#x.connect("calinput", "noelsets", r"\text{inpc,..,ier}")
#x.connect("calinput", "sectionprints", r"\text{inpc,..,ier}")
#x.connect("calinput", "materials", r"\text{inpc,..,ier}")
#x.connect("calinput", "nodes", r"\text{inpc,..,ier}")
#x.connect("calinput", "nodeprints", r"\text{inpc,..,ier}")
#x.connect("calinput", "objective", r"\text{inpc,..,flag}")
#x.connect("calinput", "outputs", r"\text{inpc,..,ier}")
#x.connect("calinput", "solidsections", r"\text{inpc,..,ier}")
#x.connect("calinput", "statics", r"\text{inpc,..,ier}")
#x.connect("calinput", "steps", r"\text{inpc,..,namtot}")
#x.connect("calinput", "surfaces", r"\text{inpc,..,ier}")
#x.connect("calinput", "gen3delem", r"\text{kin,..,jobname}")

# sensitivity.c --------------------------------#
x.connect("sensitivity", "mafillsmmainse", r"\text{co,..,eleVol}")
x.connect("mafillsmmainse","mafillsmse", r"\text{co1,..,eleVol1}")
x.connect("mafillsmse","ec3dse", r"\text{co1,..,zcg}")


#linstatic.c ------------------------------------#
x.connect("linstatic","mafillsmmain", r"\text{co,..,penal}")
x.connect("linstatic","mafillsmas", r"\text{co,..,network}")
x.connect("linstatic","spooles", r"\text{ad,..,nzs}")

# filterVector.c -------------------------------#
x.connect("mafillsmmainVectorfilter", "elementcpuload", r"\text{neapar,..,cpus}")
x.connect("filterVector", "mafillsmmainVectorfilter", r"\text{ipkon,..,q}")
x.connect("mafillsmmainVectorfilter", "mafillsmmainVectorfiltermt", r"\text{ipkon,..,q}")
x.connect("mafillsmmainVectorfiltermt", "mafillsmvectorfilter", r"\text{ne1,..,q1}")

# First tier connections ---------------------#
x.connect("main", "readinput", r"\text{jobname,..,nuel}")
x.connect("main", "allocation", r"\text{nload,..,nef}")
x.connect("main", "calinput", r"\text{cp,..,veloo}")
x.connect("main", "rho", r"\text{design,ne}")
x.connect("main", "writeheading", r"\text{jobnamec,..,nheading}")
x.connect("main", "tiedcontact", r"\text{ntie,..,kind2}")
x.connect("main", "cascade", r"\text{ipompc,..,ithermal}")
x.connect("main", "mastruct", r"\text{nk,..,network}")
x.connect("main", "densityfilter", r"\text{co,..,fnnzassumed}")
x.connect("main", "filterVector", r"\text{ipkon,..,qfilter}")
x.connect("main", "linstatic", r"\text{co,...,pstiff}")
x.connect("main", "nonlingeo", r"\text{co,..,pstiff}")
x.connect("main", "sensitivity", r"\text{co,..,eleVol}")




# Output -------------------------------------#
x.add_output("ec3dse", r"\nabla *", side=LEFT)
#---------------------------------------------#


#x.connect("adsb", "openap", r"\text{Mission parameters}")

#x.connect("CIRIUM", "engine", r"\text{Engine variant}")

#x.connect("openap", "aero", r"\text{Mission parameters, mass fraction}")
#x.connect("openap", "kin", r"\text{Mission parameters}")
#x.connect("openap", "engine", r"\text{Mission parameters, mass fraction}")

#x.connect("aero", "traj", r"\text{Lift, Drag}")
#x.connect("kin", "traj", r"\text{R.O.C, T.O. speed... }")
#x.connect("engine", "traj", r"\text{Avaiable thrust}")


#x.add_output("traj", r"\text{Fuel burn}", side=LEFT)

x.write("calTop")
