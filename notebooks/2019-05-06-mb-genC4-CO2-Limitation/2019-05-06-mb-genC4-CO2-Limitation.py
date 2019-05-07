
# coding: utf-8

# # Effect of CO2 limitation on the C4 mode
# 

# ## 0. Initialization

# In[1]:


#Import sys
import sys 
sys.path.append("../src/") 

#Import init for initialisation & loading user-defined functions
from init_fba import *

#Initialize notebook settings
theNotebook = '2019-05-06-mb-genC4-CO2-Limitation'
init_notebook(theNotebook)

#load sbml model
c3_model = load_sbml_model()


# ## 1. C3 Model

# ### 1.1 Constraints

# In[2]:


#CONSTRAINT: Set flux of all export reaction to zero
for r_obj in c3_model.reactions:
    r_id = r_obj.id
    if r_id[0:2] == "Ex":
        r_obj.bounds = (0.,0.)

#CONSTRAINT: Divergent fluxes of export and import reactions
set_bounds('Im_CO2', (-inf, inf), c3_model)
set_bounds('Im_H2O', (-inf, inf), c3_model)
set_bounds('Im_H2S', (0.,0.), c3_model)
set_bounds('Im_NH4', (0., 0.), c3_model)
set_bounds('Im_NO3', (0., inf), c3_model)
set_bounds('Im_Pi', (0., inf), c3_model)
set_bounds('Im_SO4', (0., inf), c3_model)
set_bounds('Ex_O2', (-inf, inf), c3_model)
set_bounds('Ex_Suc', (0., inf), c3_model)
set_bounds('Ex_starch', (0., inf), c3_model)
set_bounds('Ex_AA', (0., inf), c3_model)

#CONSTRAINT: 
set_bounds('G6PDH_h', (0.,0.), c3_model)
set_bounds('PPIF6PK_c', (0,0.), c3_model)

#CONSTRAINT: max. photon consumption auf C4 plants 1000 μE
set_bounds('Im_hnu', (0, 1000), c3_model)

#CONSTRAINT: CO2 uptake rate in C4 plants is higher, about 40 μmol/(m2*s)
f_CO2 = 40 #[μmol/(m2*s)] 
set_bounds('Im_CO2', (0, f_CO2), c3_model)


# In[3]:


#CONSTRAINT: Maintenace cost

atp_cost_L3_m = 0.009111187245501572 #Mitochondria-L3-ATP Cost [µmol*s-1*m-2]
atp_cost_L3_h = 0.15270708327974447 #Chloroplast-L3-ATP Cost [µmol*s-1*m-2]
atp_cost_L3_p = 0.0076669066992201855 #Peroxisome-L3-ATP Cost [µmol*s-1*m-2]
atp_cost_L3_c = 0.042683072918274702 #Cytosl/Other-L3-ATP Cost [µmol*s-1*m-2]

set_fixed_flux('NGAM_c',atp_cost_L3_c + atp_cost_L3_p, c3_model)
set_fixed_flux('NGAM_m',atp_cost_L3_m, c3_model)
set_fixed_flux('NGAM_h',atp_cost_L3_h, c3_model)


# In[4]:


#CONSTRAINT: fluxes through the chloroplastic NADPH dehydrogenase and plastoquinol oxidase were set to zero 
#because the contributions of NADPH dehydrogenase (Yamamoto et al., 2011) and plastoquinol oxidase 
#(Josse et al., 2000) to photosynthesis are thought to be minor.
set_bounds('AOX4_h',(0,0), c3_model)
set_bounds('iCitDHNADP_h',(0,0), c3_model)


# In[5]:


#CONSTRAINT: NTT is only active at night
set_fixed_flux('Tr_NTT',0, c3_model)


# In[6]:


#CONSTRAINT: No uncoupled pyruvate transport
set_bounds('Tr_Pyr1',(0,0), c3_model)
set_bounds('Tr_Pyr2',(0,0), c3_model)


# ## 2. C4 Model

# ### 2.1 Compose Model

# #### 2.1.1 Two copies of genC3 model and exchange reactions

# In[7]:


c4_model = cobra.Model('c4_model')

cell_types = ['M', 'B']

#duplicate metabolites
for m in c3_model.metabolites:
    for cell in cell_types:
        m_dt = cobra.Metabolite('['+cell+']_'+m.id, name = m.formula, compartment = m.compartment)
        c4_model.add_metabolites([m_dt])

#duplicate reactions
for r_c3_obj in c3_model.reactions:
    for cell in cell_types:
        r_c4_obj = cobra.Reaction('['+cell+']_'+r_c3_obj.id)
        r_c4_obj.name = r_c3_obj.name
        r_c4_obj.subsystem = r_c3_obj.subsystem
        r_c4_obj.bounds = r_c3_obj.bounds
        c4_model.add_reaction(r_c4_obj)
        r_c4_obj.add_metabolites({'['+cell+']_'+m_c3_obj.id: r_c3_obj.get_coefficient(m_c3_obj) for m_c3_obj in r_c3_obj.metabolites})
        
        
#metabolites excluded from M/BS exchange
no_transport = ['NO3','NO2', 'O2','Na', 'H2S', 'SO4',
                'H2O','FBP','F26BP','DPGA','H','ACD','AC','M_DASH_THF', '5M_DASH_THF', 'H_DASH_Cys', 'aH_DASH_Cys', 'ORO', 'DHO',
                'GABA','A_DASH_Ser','PRPP','AD','THF','DHF','ADN','Mas','CoA','GluP',
                'A_DASH_CoA','cellulose1','cellulose2','cellulose3','starch1',
                'starch2','starch3','TRXox','TRXrd','Glu_DASH_SeA','T6P','aMet',
                'PPi', 'P5C', 'NH4', 'Pi', 'CO2', 'OAA','HCO3', 
                'UTP', 'UDP', 'UDPG', 'ATP', 'ADP', 'AMP', 'IMP', 'XMP', 
                'GTP', 'GDP', 'GMP', 'OMP', 'UMP', 'CTP', 'GDP', 'CDP', 'dADP', 
                'dCDP', 'dGDP', 'dUDP', 'dUTP', 'dUMP', 'dTMP', 'dTDP', 'GTP', 
                'dATP', 'dCTP', 'dGTP', 'dTTP', 'NAD', 'NADH', 'NADP', 'NADPH']

#add M/BS exchange reactions
L_r_transport = []
for m_c3_obj in c3_model.metabolites:
    if m_c3_obj.id[-1:] == 'c' and m_c3_obj.id[:-2] not in no_transport:
        r_c4_obj = cobra.Reaction('[MB]_'+m_c3_obj.id)
        r_c4_obj.name = '[MB]_'+m_c3_obj.id
        r_c4_obj.subsystem = 'Exchange'
        r_c4_obj.bounds = (-inf, inf)
        c4_model.add_reaction(r_c4_obj)
        r_c4_obj.add_metabolites({'[M]_'+m_c3_obj.id: -1,'[B]_'+m_c3_obj.id: 1 })
        L_r_transport.append('[MB]_'+m_c3_obj.id)


# #### 2.1.2 Adaptations for second unconstrained rubisco population in the bundle sheath

# In[8]:


#CONSTRAINT: Add external CO2 species in bundle sheat
#(the original CO2 species is treated as internal CO2)
m_list_CO_Ex= ['[B]_CO2_ex_c','[B]_CO2_ex_h','[B]_CO2_ex_m']

for m_id in m_list_CO_Ex:
    m_obj = cobra.Metabolite(m_id)
    c4_model.add_metabolites(m_obj)

#CONSTRAINT: Copy all reactions using internal CO2 and exchange internal with external CO2 in the copied recactions
r_list_CO_Ex = ['Tr_CO2h', 'RBC_h']

for r_id in r_list_CO_Ex:
    r_obj = c4_model.reactions.get_by_id('[B]_'+r_id)
    r_obj_Ex = cobra.Reaction(r_obj.id+'_Ex')
    r_obj_Ex.name = r_obj.id+'_Ex'
    r_obj_Ex.subsystem = r_obj.subsystem
    r_obj_Ex.bounds = r_obj.bounds
    c4_model.add_reaction(r_obj_Ex)
    r_obj_Ex.add_metabolites({m_obj.id if not m_obj.id[:-2] == '[B]_CO2' else '[B]_CO2_ex'+m_obj.id[-2:]: r_obj.get_coefficient(m_obj) 
                                  for m_obj in r_obj.metabolites})

#CONSTRAINT: CO2 exchange between mesophyll and bundle sheat
r_c4_obj = cobra.Reaction('[MB]_CO2_c')
r_c4_obj.name = '[MB]_CO2_c'
r_c4_obj.subsystem = 'Exchange'
r_c4_obj.bounds = (-inf, inf)
c4_model.add_reaction(r_c4_obj)
r_c4_obj.add_metabolites({'[M]_CO2_c': -1,'[B]_CO2_ex_c': 1 })
L_r_transport.append('[MB]_CO2_c')


# ### 2.3 Constraints

# In[9]:


#CONSTRAINT: No CO2 uptake in bundle sheat cells due to suberin layer in cell membranes
set_fixed_flux('[B]_Im_CO2',0, c4_model)


# In[10]:


#CONSTRAINT: Output of sucrose : total amino acid and sucrose : starch
set_fixed_flux_ratio({'[B]_Ex_Suc':2.2,'[B]_Ex_AA':1.0}, c4_model)
set_fixed_flux_ratio({'[B]_Ex_Suc':1.0,'[B]_Ex_starch':1.0}, c4_model)

set_fixed_flux_ratio({'[M]_Ex_Suc':2.2,'[M]_Ex_AA':1.0}, c4_model)
set_fixed_flux_ratio({'[M]_Ex_Suc':1.0,'[M]_Ex_starch':1.0}, c4_model)


# In[11]:


#Reaction variables for light uptake
B_Im_hnu = c4_model.reactions.get_by_id("[B]_Im_hnu")
M_Im_hnu = c4_model.reactions.get_by_id("[M]_Im_hnu")

#CONSTRAINT: Total Photon uptake limited to 1000 µE
const_hnu_sum = c4_model.problem.Constraint( B_Im_hnu.flux_expression + M_Im_hnu.flux_expression,
                                        lb = 0, ub = 1000)

c4_model.add_cons_vars(const_hnu_sum)

#CONSTRAINT: Total Photon uptake by bundle sheath must be less equal than in mesophyll
const_hnu_ratio = c4_model.problem.Constraint( M_Im_hnu.flux_expression - B_Im_hnu.flux_expression,
                                        lb = 0, ub = 1000)

c4_model.add_cons_vars(const_hnu_ratio)


# In[12]:


#CONSTRAINT: oxygenation : carboxylation = 1 : 3
set_fixed_flux_ratio({'[B]_RBC_h_Ex':3,'[B]_RBO_h':1}, c4_model)
set_fixed_flux_ratio({'[M]_RBC_h':3,'[M]_RBO_h':1}, c4_model)


# ## 3 FBA

# In[13]:


#Dictionary defining the experiments according to the CO2-limitation
D_exp = {value: value for value in np.arange(0,41,5)}

#Reaction Variables
B_Ex_Suc = c4_model.reactions.get_by_id("[B]_Ex_Suc")
B_RBO = c4_model.reactions.get_by_id("[B]_RBO_h")
M_RBO = c4_model.reactions.get_by_id("[M]_RBO_h")
M_Im_CO2 = c4_model.reactions.get_by_id("[M]_Im_CO2")

#Set FBA solver
c4_model.solver = "glpk"

#Initialize dictionary to store results
D_fba={}

#Run every FBA experiment

for value in sorted(D_exp.keys()): #iterate over proportions of carboxylation
    
    #Add CO2-limitation constraint 
    const_CO2 = c4_model.problem.Constraint( M_Im_CO2.flux_expression,
                                        lb = value, ub = value)
    c4_model.add_cons_vars(const_CO2)
    
    #Optimize/Maximize sucrose output
    B_Ex_Suc.objective_coefficient = 1.
    c4_model_copy = c4_model.copy()
    result_fba = c4_model_copy.optimize('maximize')
    del c4_model_copy
    c4_model.objective =[]
    set_fixed_flux(B_Ex_Suc.id,result_fba.fluxes[B_Ex_Suc.id], c4_model)
    
    #Optimize/Minimize oxygenation rate by rubisco (set True)
    if True:
            B_RBO.objective_coefficient = 1.
            M_RBO.objective_coefficient = 1.
            c4_model_copy = c4_model.copy()
            result_fba = c4_model_copy.optimize('minimize')
            del c4_model_copy
            c4_model.objective =[]
            set_fixed_flux(B_RBO.id,result_fba.fluxes[B_RBO.id],c4_model)            
            set_fixed_flux(M_RBO.id,result_fba.fluxes[M_RBO.id],c4_model)
    
        
    #Optimize/Minimize total flux
    if result_fba.status == 'optimal': 
        c4_model_copy = c4_model.copy()
        result_pfba = cobra.flux_analysis.parsimonious.pfba(c4_model_copy)
        D_fba[value] = result_pfba.fluxes
        del c4_model_copy 
    
    #Reset reaction bounds
    set_bounds(B_RBO.id,(0,inf),c4_model)
    set_bounds(M_RBO.id,(0,inf),c4_model)
    set_bounds(B_Ex_Suc.id,(0,inf),c4_model)

    #Remove CO2-limitation constraint
    c4_model.remove_cons_vars(const_CO2)


# ## 4 Figures

# In[14]:


xaxis_title = 'CO2 Uptake [µmol/s/m2]'
save_fig = False


# In[15]:


create_bar_plot_rxn(D_fba, D_exp,
                    {'[M]_Im_hnu':'Mesophyll','[B]_Im_hnu': 'Bundle sheath'},
                    'Photon Uptake', xaxis_title,
                    stacked = True, save_fig=save_fig)


# In[16]:


create_bar_plot_met(D_fba, D_exp,
                    {'[M]_PSII_h':'Mesophyll','[B]_PSII_h': 'Bundle sheath'},'hnu_h',
                    'Photon Uptake by PSII', xaxis_title, c3_model,
                    stacked = True, save_fig=save_fig)


# In[17]:


create_bar_plot_met(D_fba, D_exp,
                    {'[M]_PSI_h':'Mesophyll','[B]_PSI_h': 'Bundle sheath'},'hnu_h',
                    'Photon Uptake by PSI', xaxis_title, c3_model,
                    stacked = True, save_fig=save_fig)


# In[18]:


create_bar_plot_met(D_fba, D_exp,
                    {'[M]_cplx5_m':'Mesophyll','[B]_cplx5_m': 'Bundle sheath'},'ATP_m',
                    'ATP Synthesis Mitochondria', xaxis_title, c3_model,
                    stacked = True, save_fig=save_fig)


# In[19]:


create_bar_plot_met(D_fba, D_exp,
                    {'[M]_ATPase_h':'Mesophyll','[B]_ATPase_h': 'Bundle sheath'},'ATP_h',
                    'ATP Synthesis Chloroplast', xaxis_title, c3_model,
                    stacked = True, save_fig=save_fig)


# In[20]:


create_bar_plot_rxn(D_fba, D_exp,
                    {'[M]_RBC_h':'Mesophyll (constrained)','[B]_RBC_h_Ex':'Bundle sheat (constrained)',
                 '[B]_RBC_h': 'Bundle sheat (unconstrained)'},
                    'Carboxylation by Rubisco', xaxis_title,
                    stacked = True, save_fig=save_fig)


# In[21]:


create_bar_plot_rxn(D_fba, D_exp,
                    {'[B]_MalDH4_h': 'NADP-ME', '[B]_MalDH2_m': 'NAD-ME', '[B]_PEPC1_c': 'PEPCK',},
                    'Decarboxylation Enzymes', xaxis_title,
                    c=True, stacked = True, save_fig=save_fig)


# In[22]:


create_bar_plot_rxn(D_fba, D_exp,
                    {'[M]_PEPC2_c': 'Mesophyll', '[B]_PEPC2_c': 'Bundle sheath'},
                    'PEPC', xaxis_title,
                     stacked = True, save_fig=save_fig)


# In[23]:


create_bar_plot_rxn(D_fba, D_exp,
                    {'[M]_PyrPiDK_h': 'Mesophyll', '[B]_PyrPiDK_h': 'Bundle sheath'},
                    'PyrPiDK', xaxis_title,
                     stacked = True, save_fig=save_fig)


# In[24]:


plot_transport(D_fba, D_exp, L_r_transport, xaxis_title, save_fig=save_fig)

