# Reaction fluxes and compound usage
## Introduction
#In this notebook we will highlight features in _PyFBA_ that allows the user to obtain a model's reaction distribution after performing flux-balance analysis.
#**Data used in this notebook**
##  _Citrobacter sedlakii_ genome-scale metabolic model
##  LB media
import sys
import os

# Add local PyFBA to sys path so we can use the developer version
sys.path.insert(0, os.path.expanduser("~") + "/PyFBA/")
import PyFBA

## Load metabolic model

#_PyFBA_ contains two methods to load and save metabolic models.
#1. Metabolic model directories
#2. SBML files

#For this example, we have the _Citrobacter sedlakii_ model stored in a model directory.
# Specify location of metabolic model
modelDir = "citrobacter_sedlakii/"
orgName = "citrobacter_sedlakii"

model = PyFBA.model.load_model(modelDir, orgName)
print("Number of reactions", len(model.reactions))
print("Number of compounds", len(model.compounds))

# Test its ability to grow in the LB media "ArgonneLB.txt"
status, biomassFlux, growth = model.run_fba(media_file="MOPS_NoC_D-Glucose.txt")
print("Biomass flux: {:.3f}".format(biomassFlux))
print("Growth:", growth)

## Obtain flux distribution from reactions

#_PyFBA_ provides the capability to obtain the flux distribution resulting from a flux-balance analysis process. Using the `PyFBA.model.model_reaction_fluxes()` function, each reaction and their flux is returned as a Dictionary where:
#- **KEY:** ModelSEED Reaction ID (String)
#- **VALUE:** Reaction flux (float)

# Get number of reactions with a non-zero flux
nonzeroRxns = {r: flux for r, flux in rxnFluxes.items() if flux != 0}
print("Number of non-zero flux reactions:", len(nonzeroRxns), "out of", len(rxnFluxes))
print("Percent of reactions with non-zero flux: {:.2f}%".format(len(nonzeroRxns) / len(rxnFluxes) * 100))


### Uptake and secretion reactions

#The uptake and secretion reactions are designated as "reaction sinks". In practice they are indicated as **EX** reactions or with the SBML attribute **boundaryCondition** set to _True_. These reaction sinks are necessary to satisfy the steady-state assumption of flux-balance analysis.

#In _PyFBA_, we add these reactions when flux-balance analysis is executed. These are added based on all external compounds within the system:
#1. Media compounds
#2. External compounds from transports

#The flux values for these reaction sinks indicate the compounds transported into the cell from the media (**negative flux for uptake**) and the compounds secreted to the environment (**positive values for secretion**).


# Get uptake/secretion reactions
usRxns = {r: flux for r, flux in rxnFluxes.items() if "UPTAKE" in r}
print("Number of uptake/secretion reactions:", len(usRxns))

## Compound usage counts

#In addition to the reaction fluxes, _PyFBA_ provides the capability to obtain the usage of each compound within the metabolic model. This is calculated by taking the sum of all reaction fluxes that either consume or produce the compound. It is important to note that flux-balance analysis assumes a steady-state where the change in individual compound concentrations is zero. However, what these compound counts will indicate is their overall "usage" to obtain the optimal biomass production.

#Using the `PyFBA.model.compute_compound_counds()` function, a multi-level Dictionary is returned.

#The **first level** of data are the model components location, for example:
#1. "e" for extracelluar
#2. "c" for cytoplasm

#---
#The **second level** of data are the PyFBA Compound String representations, for example:
#- `Glycine (location: e)`

#---
#The **third level** of data is a Dictionary containing two keys:
#1. "count" - The value for count is a `float` for the final sum of all reaction fluxes, for example:
#   - `{"Glycine (location: e)" : {"count": 1000.0}, "reactions": {...}}`
#2. "reactions" - The value for reactions is a Dictionary discussed below.

#---
#The **fourth level** of data is the "reactions" Dictionary that contains each reaction consuming or producing the compound and their flux values. These are the flux values used to calculate the "count" in the level above. For example:
#- `"reactions": {"rxn0001": 1000.0, "rxn0002": 0.0, "rxn0003": -1000.0, ...}`

#---
#One important note regarding Uptake/Secretion reactions is that those fluxes are not included in the summation because these reactions are artificial. What this results in is non-zero counts for some external compounds. Specifically, **negative counts** refer to external compounds being consumed (UPTAKE) and **positive counts** refer to external compounds being secreted (SECRETION).


# Print out uptake/secretion fluxes
# Negative values indicate usage by the cell (UPTAKE)
# Positive values indicate secretion by the cell (SECRETION)
for r, flux in usRxns.items():
    threshold = 0.000001
    # Skip very small numbers
    if flux > -threshold and flux < threshold:
        continue
        
    print("{}:\t{:.3f}".format(r, flux))

# The below cell runs through all of the media for our bacteria and creates two txt files labeled table_1 and table_2 followed by the media evaluated.

# Table 1 is formated as follows
# Reaction number ; the flux of that reaction ; the molecules involved in that reaction.

# Table 2 is formated as follows
# Compound (location of the compound either e or c ) ; the "count" of the compound as described in cell 10 ; the reactions that compound is involved in ; the clux of each reaction the compound was involved in.

# Finally there is a file starting with Compounds_ followed by the media grown in
# this is a text file with the compounds names one per each line.

fin = (open("./citrobacter_sedlakii/citrobacter_sedlakii.gfmedia"))
numlines = sum(1 for line in fin)
fin.close()
fin = (open("./citrobacter_sedlakii/citrobacter_sedlakii.gfmedia"))
medialist = []
for i in range(numlines):
        file = (fin.readline())
        file = file.replace("\n", '')
        medialist.append(file)
medianames = []
for i in medialist:
    newfile = i.replace(".txt","")
    medianames.append(newfile)
print(medialist)
print(medianames)
counter = 0 
for file in medialist:
    rxnFluxes = PyFBA.model.model_reaction_fluxes(model=model, media_file=file)
    cpdCounts = PyFBA.model.compute_compound_counts(model=model, media_file=file)
    fo = open("table_1" + str(medianames[counter]) + ".txt", "w")
    rxnF_list = list(rxnFluxes.keys())
    for i in rxnF_list:
        if "UPTAKE" in i:
            j = "UPTAKE"
        elif "BIOMASS" in i:
            j = "BIOMASS_EQN"
        else:
            j = str(model.reactions[i].equation) 
        fo.write( str(i) + ";" + str(rxnFluxes[i]) + ";" + j + "\n")
    fo.close()
    fo = open("table_2" + str(medianames[counter]) + ".txt", "w")
    cpds_list = []
    for cpd, data in cpdCounts["e"].items():
        if  abs(data["count"]) >= 0.0000:
            fo.write("{};{}".format(cpd, data["count"]) + "; " + str(list(cpdCounts["e"][cpd]["reactions"].keys())) + "; " + str(list(cpdCounts["e"][cpd]["reactions"].values())) + "\n")
        fo.write("\n \n")
        cpds_list.append(cpd)
    for cpd, data in cpdCounts["c"].items():
        if abs(data["count"]) >= 0.0000:
            fo.write("{};{}".format(cpd, data["count"]) + "; " + str(list(cpdCounts["c"][cpd]["reactions"].keys())) + "; " + str(list(cpdCounts["c"][cpd]["reactions"].values())) + "\n")
        fo.write("\n \n")
        cpds_list.append(cpd)
    fo.close()
    fo = open('compounds_' +  str(medianames[counter]) + ".txt",'w')
    cpds_list = list(set(cpds_list))
    for cpd in cpds_list:
        fo.write(cpd + '\n')
    fo.close()
    counter = counter + 1
