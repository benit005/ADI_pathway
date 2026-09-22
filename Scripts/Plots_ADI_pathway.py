#!/usr/bin/env python

"""
Author Sara Benito. 
Script to generate the plots for Fig 5S (A/B) of the publication
Cobrapy and functions are used from https://github.com/opencobra/m_model_collection/
Python 3.14 is used as the programming language with the cobrapy version 0.32.1
"""

# Import statements

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.lines import Line2D


if __name__ == "__main__":
	file = "../Output/Simulations_ADI_pathway.xlsx" #Output file with modelling results
	conditions = ["Glc_WT","Glc_Arg_WT","Glc_ADI","Glc_Arg_ADI"] # Use for stacked bar plots and FVA ranges
	products = ["EX_glu__L_e","EX_pyr_e","EX_ac_e","EX_lac__D_e","EX_for_e"] #Products shown on Fig5SB
	products_stacked = ["EX_orn_e","EX_citr__L_e"] #Products shown on the stacked bar plot Fig5SA
	product_names = {"EX_orn_e": "Ornithine","EX_citr__L_e": "Citrulline","EX_glu__L_e": "Glutamate","EX_pyr_e": "Pyruvate",\
	"EX_ac_e": "Acetate","EX_lac__D_e": "Lactate","EX_for_e": "Formate"} #Dictionary with all the products and their names
	
	#Read excel file 
	# Sheet 1: FBA fluxes
	df_flux = pd.read_excel(file,sheet_name=0,index_col=0)
	# Sheet 3: FVA
	df_fva = pd.read_excel(file,sheet_name="ADI_pathway_FVA",index_col=0)
	
	#Prepare data for Fig5SA
	biomass_reaction = ("('BIOMASS_Ec_iJO1366_WT_53p95M', ""'E. coli biomass objective function (iJO1366) - WT - with 53.95 GAM estimate')")
	Glucose_uptake=10
	growth_rates = df_flux.loc[biomass_reaction,conditions]

	# Calculate biomass yield to show on the same plot
	biomass_yield = growth_rates / Glucose_uptake * 1000

	fva_data = []
	for product in products + products_stacked:
		for condition in conditions:
			fva_data.append({"condition": condition,"product": product_names[product],\
			"minimum": df_fva.loc[product,f"{condition}_min"],"maximum": df_fva.loc[product,f"{condition}_max"]})

	fva_df = pd.DataFrame(fva_data)
	#Dataframe for the stacked plot
	
	stacked_df = fva_df[fva_df["product"].isin([product_names[p] for p in products_stacked])].pivot(
	index="condition",
	columns="product",
	values="maximum")
	stacked_df = stacked_df[[ "Ornithine","Citrulline"]]
	stacked_df = stacked_df.loc[conditions]

	# Define the colors for the products shown on the stacked bar plot 
	product_colors = {"Ornithine": "cadetblue","Citrulline": "mediumaquamarine"}

	# Color for biomass yield marker
	biomass_color = "#FBDD7E"
	#sns.set_theme(style="whitegrid")
	
	#Plot stacked bar plots with biomass yield
	
	fig, ax1 = plt.subplots(figsize=(6, 6))

	bottom = np.zeros(len(stacked_df))

	for product in stacked_df.columns:
		ax1.bar(stacked_df.index,stacked_df[product],bottom=bottom,label=product,color=product_colors[product])
		bottom += stacked_df[product].values

	ax1.set_ylim(0, 25)
	ax1.set_ylabel("Maximum production rate (mmol gDW$^{-1}$ h$^{-1}$)")
	
	#Include the biomass yield on the secondary y axis
	ax2 = ax1.twinx()
	x = range(len(conditions))

	ax2.plot(x,biomass_yield.values,marker="D",linestyle="None",markersize=9,markeredgecolor="black",markeredgewidth=1.5,\
	color=biomass_color,label="Biomass yield")
	ax2.set_ylabel("Biomass yield (g mol$^{-1}$ glc)")
	# Legends
	ax1.legend(bbox_to_anchor=(0.03, 1),loc="upper left")
	ax2.legend(bbox_to_anchor=(0.03, 0.9),loc="upper left")

	plt.tight_layout()
	#Save figure 
	plt.savefig("../Output/Fig5SA.png", dpi=300)
	plt.savefig("../Output/Fig5SA.svg", dpi=600)
	#plt.show()

	#Prepare Fig5SB - FVA results for the selected products

	fig, ax = plt.subplots(figsize=(8, 6))
	product_list = [product_names[p]for p in products]
	n_conditions = len(conditions)
	offsets = [-0.15, -0.05, 0.05, 0.15]
	#Colors
	condition_colors = {"Glc_WT": "#8C000F", "Glc_Arg_WT": "#DC143C","Glc_ADI": "#008080","Glc_Arg_ADI": "#06C2AC"}

	for i, condition in enumerate(conditions):
		condition_data = fva_df[fva_df["condition"] == condition]

		# Color assigned to this condition
		color = condition_colors[condition]

		for j, product in enumerate(product_list):
			row = condition_data[condition_data["product"] == product].iloc[0]
			y = j + offsets[i]
			ax.plot([row["minimum"], row["maximum"]],[y, y],color=color)
			ax.plot(row["minimum"],y,marker="o",color=color,markersize=5)
			ax.plot(row["maximum"],y,marker="o",color=color,markersize=5)


	ax.set_yticks(range(len(product_list)))
	ax.set_yticklabels(product_list)

	ax.set_xlabel("Minimum and maximum production rate (mmol gDW$^{-1}$ h$^{-1}$)")
	ax.set_ylabel("Products")

	legend_handles = [Line2D([0],[0],color=condition_colors[condition],linewidth=3,marker="o",markersize=5,label=condition)
		for i, condition in enumerate(conditions)]
	ax.legend(
		handles=legend_handles,
		title="Condition",
		bbox_to_anchor=(1.02, 1),
		loc="upper left")
	#Save the plot
	plt.tight_layout()
	plt.savefig("../Output/Fig5SB.png", dpi=300)
	plt.savefig("../Output/Fig5SB.svg", dpi=600)
	#plt.show()
