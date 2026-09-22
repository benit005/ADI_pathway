#!/usr/bin/env python

"""
Author Sara Benito. 
Script to perform FBA/FVA simulations to assess E. coli's growth with and without addition of  
the ADI pathway on different carbon sources.  
Cobrapy and functions are used from https://github.com/opencobra/m_model_collection/
Python 3.14 is used as the programming language with the cobrapy version 0.32.1
"""
# Import statements

import cobra
from cobra import Model, Reaction, Metabolite
import pandas
from pandas import DataFrame
import numpy as np
import pandas as pd
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)
from cobra.flux_analysis import flux_variability_analysis
import statistics as stats
from statistics import stdev
import cobra.util.solver as sutil
from cobra.flux_analysis import pfba


def add_ADI_pathway(Model):
	"""
	Function to add add the reactions to
	complement the ADI pathway in E.coli's model

	Input:
	model: model, cobrapy model structure

	Outpu:
	model, cobrapy model structure with the added pathway
	"""

# Define the reactions

	r1 = Reaction("argA")
	r1.name = "Arginine deiminase"
	r1.lower_bound = 0.  
	r1.upper_bound = 1000 
	r1.add_metabolites({
	model.metabolites.get_by_id("arg__L_c"): -1.0,
	model.metabolites.get_by_id("h2o_c"):-1.0,
	model.metabolites.get_by_id("nh4_c"): +1.0, 
	model.metabolites.get_by_id("citr__L_c"): +1.0})
	
	r2 = Reaction("TCITRe")
	r2.name = "Citrulline transport"
	r2.lower_bound = 0.  
	r2.upper_bound = 1000 
	r2.add_metabolites({
	model.metabolites.get_by_id("citr__L_c"): -1.0,
	model.metabolites.get_by_id("h_c"): -1.0,
	model.metabolites.get_by_id("h_e"): 1.0,
	model.metabolites.get_by_id("citr__L_e"): +1.0})	

	r3 = Reaction("EX_citr__L_e")
	r3.name = "Citrulline exchange"
	r3.lower_bound = 0.  # This is the default
	r3.upper_bound = 1000. # This is the default
	r3.add_metabolites({
	model.metabolites.get_by_id("citr__L_e"): -1.0})	

	model.add_reactions({r1,r2,r3})	


	return model

def no_ADI_pathway(Model):
	"""
	Function to add the reactions of the ADI pathway
	in the model (inactive) to have the same list of reactions 
	in the different conditions

	Input:
	model: model, cobrapy model structure

	Outpu:
	model, cobrapy model structure with the added pathway
	"""

# Define the reactions

	r1 = Reaction("argA")
	r1.name = "Arginine deiminase"
	r1.lower_bound = 0.  
	r1.upper_bound = 0. # Conditions without the ADI pathway
	r1.add_metabolites({
	model.metabolites.get_by_id("arg__L_c"): -1.0,
	model.metabolites.get_by_id("h2o_c"):-1.0,
	model.metabolites.get_by_id("nh4_c"): +1.0, 
	model.metabolites.get_by_id("citr__L_c"): +1.0})
	
	r2 = Reaction("TCITRe")
	r2.name = "Citrulline transport"
	r2.lower_bound = 0.  # 
	r2.upper_bound = 0. # No ADI pathway
	r2.add_metabolites({
	model.metabolites.get_by_id("citr__L_c"): -1.0,
	model.metabolites.get_by_id("h_c"): -1.0,
	model.metabolites.get_by_id("h_e"): 1.0,
	model.metabolites.get_by_id("citr__L_e"): +1.0})	

	r3 = Reaction("EX_citr__L_e")
	r3.name = "Citrulline exchange"
	r3.lower_bound = 0.  # 
	r3.upper_bound = 0. # No ADI pathway
	r3.add_metabolites({
	model.metabolites.get_by_id("citr__L_e"): -1.0})	

	model.add_reactions({r1,r2,r3})

	return model


if __name__ == "__main__":
	model=cobra.io.read_sbml_model("../Data/iJO1366.xml") #E.coli's model
	solution_list=[] #List to store the FBA solutions simulated per condition
	model.reactions.EX_o2_e.bounds=[0,0] #Simulate anaerobic conditions
	model.reactions.ATPM.bounds=[8.4,1000] #Set a default lb for ATPM
	
	products = ["EX_orn_e","EX_citr__L_e","EX_glu__L_e","EX_pyr_e","EX_ac_e", "EX_lac__D_e","EX_for_e"] #List of products to track

	# 1 Assess growth on glucose with/o arginine and with/o ADI pathway
	
	with model:
		#Glucose as carbon source (no ADI)

		no_ADI_pathway(model) #We just add the Adi reactions with no flux to have the same list of reactions
		model.reactions.EX_glc__D_e.bounds=[-10,-10] #
		solution_glc = model.optimize() #fba
		pfba_sol_glc = cobra.flux_analysis.pfba(model) #pFBA
		fva_glc = flux_variability_analysis(model,reaction_list=products,fraction_of_optimum=0.9) #FVA

		#Glucose and arginine as carbon source (no ADI)
		
		model.reactions.EX_arg__L_e.bounds=[-20,-20] #To keep the experimental glc/arg ratio
		solution_glc_arg = model.optimize()
		pfba_sol_glc_arg = cobra.flux_analysis.pfba(model)
		fva_glc_arg = flux_variability_analysis(model,reaction_list=products,fraction_of_optimum=0.9)

	with model:	
		#Glucose as carbon source with ADI pathway
		add_ADI_pathway(model)
		model.reactions.EX_glc__D_e.bounds=[-10,-10]
		solution_glc_adi = model.optimize()
		pfba_sol_glc_adi = cobra.flux_analysis.pfba(model)	
		fva_glc_adi = flux_variability_analysis(model,reaction_list=products,fraction_of_optimum=0.9)
		
		#Glucose and arginine as carbon source with ADI pathway
		
		model.reactions.EX_arg__L_e.bounds=[-20,-20]
		solution_glc_arg_adi = model.optimize()
		pfba_sol_glc_arg_adi = cobra.flux_analysis.pfba(model)
		fva_glc_arg_adi = flux_variability_analysis(model,reaction_list=products,fraction_of_optimum=0.9)

					
	# 2 Assess growth on acetate and arginine with ADI pathway (evolved strain)
	with model:
		model.reactions.EX_glc__D_e.bounds=[0,0]
		model.reactions.EX_ac_e.bounds=[-20,-5]
		model.reactions.EX_arg__L_e.bounds=[-20,-5]
		model.reactions.ATPS4rpp.bounds=[0,0] #Knock out ATP synthase reaction
		#Activate ADI pathway
		add_ADI_pathway(model)
		solution_ace_arg_adi = model.optimize()
		pfba_sol_ace_arg_adi = cobra.flux_analysis.pfba(model)	
		
		#Add Dld, sdhA, and nuoG knock outs as the evolved strain
		model.genes.b2133.knock_out()
		#model.genes.b1380.knock_out() #Also involved in LDH_D (isoenzyme)
		model.genes.b2283.knock_out()
		model.genes.b0723.knock_out()
		solution_ace_arg_adi_evol = model.optimize()
		pfba_sol_ace_arg_adi_evol = cobra.flux_analysis.pfba(model)
		#print(model.summary())
		
	# Store the results and create the output file
		reactions=[]
		for i in model.reactions:
			reactions.append((i.id, i.name))
	y=[]   #Store the data form FBA
	y1=[]	#Store the data from pFBA solutions
	fva_data=[] #Store the fva results
	for reaction_id, reaction_name in reactions:
		y.append((solution_glc.fluxes[reaction_id],solution_glc_arg.fluxes[reaction_id],solution_glc_adi.fluxes[reaction_id],\
		solution_glc_arg_adi.fluxes[reaction_id],solution_ace_arg_adi.fluxes[reaction_id],solution_ace_arg_adi_evol.fluxes[reaction_id]))
		y1.append((pfba_sol_glc.fluxes[reaction_id],pfba_sol_glc_arg.fluxes[reaction_id],pfba_sol_glc_adi.fluxes[reaction_id],\
		pfba_sol_glc_arg_adi.fluxes[reaction_id],pfba_sol_ace_arg_adi.fluxes[reaction_id],pfba_sol_ace_arg_adi_evol.fluxes[reaction_id]))
	for reaction_id in products:
		fva_data.append((fva_glc.loc[reaction_id, "minimum"],fva_glc.loc[reaction_id, "maximum"],\
		fva_glc_arg.loc[reaction_id, "minimum"],fva_glc_arg.loc[reaction_id, "maximum"],\
		fva_glc_adi.loc[reaction_id, "minimum"],fva_glc_adi.loc[reaction_id, "maximum"],fva_glc_arg_adi.loc[reaction_id, "minimum"],\
		fva_glc_arg_adi.loc[reaction_id, "maximum"]))
	names=["Glc_WT","Glc_Arg_WT","Glc_ADI","Glc_Arg_ADI","Ace_Arg_ADI","Ace_Arg_ADI_evol"]
	names_fs=["Glc_WT_m","Glc_WT_std","Glc_Arg_WT_m","Glc_Arg_WT_std","Glc_ADI_m","Glc_ADI_std",\
	"Glc_Arg_ADI_m","Glc_Arg_ADI_std", "Ace_Arg_ADI_m","Ace_Arg_ADI_std","Ace_Arg_ADI_evol_m","Ace_Arg_ADI_evol_std"]
	names_fva=["Glc_WT_min","Glc_WT_max","Glc_Arg_WT_min","Glc_Arg_WT_max","Glc_ADI_min","Glc_ADI_max",\
	"Glc_Arg_ADI_min","Glc_Arg_ADI_max"]
	
	
	#Create the dataframes
	df1 = DataFrame(y, index=reactions, columns=names)
	df2=DataFrame(y1, index=reactions, columns=names)
	df3 = DataFrame(fva_data,index=products,columns=names_fva)
	writer = pd.ExcelWriter('../Output/Simulations_ADI_pathway.xlsx', engine='openpyxl')
	# #writer2=pd.ExcelWriter('Simulations_ADI4.xlsx', engine='openpyxl')
	sheet1 = 'ADI_pathway_FBA'
	sheet2 = 'ADI_pathway_pFBA'
	sheet3= 'ADI_pathway_FVA'
	#Save it in an excel file
	df1.to_excel(writer, sheet_name=sheet1)
	df2.to_excel(writer, sheet_name=sheet2)
	df3.to_excel(writer, sheet_name=sheet3)
	writer.close()

	
