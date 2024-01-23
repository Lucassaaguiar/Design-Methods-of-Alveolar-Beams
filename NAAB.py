##############################################################################
#               Programa de cálculo de vigas alveolares                      #
#                                                                            #
# Objetivo: Determinar a carga última e os modos de falhas das               #
#           vigas alveolares. Exportar resultados.                           #
# Programaddor: Lucas Alves de Aguiar                                        #
# Início: 15/09/2022                                                         #
# Revisão: 29/01/2023                                                        #  
# Atualização: 22/01/2024                                                    # 
##############################################################################

from DPAB import *
import warnings
warnings.filterwarnings('ignore')

#----------------------------------------------------------------------------#
# Main                                                                       #
#----------------------------------------------------------------------------#

# Verificação dos métodos

de         = pd.read_excel('exemples.xlsx', sheet_name='dados')         # Importando dados de planilha de Excel
n_vigas    = len(de.values[:,0])

# Método Veríssimo et. al. (2012)

A_data = np.empty((n_vigas,9), dtype = 'object_')

print('----------------------------------------')
print('      Método Veríssimo et al. (2012)    ')
print('----------------------------------------','\n')

for j in range(n_vigas):
    
    Met_A = viga_al(j)
    carga        = round(((30.72*3) * Met_A.melast * Met_A.Ix_O / (5*Met_A.dg)**3),1) 
    passo        = 0.1
    inicio       = round(((30.72/3) * Met_A.melast * Met_A.Ix_O / (20*Met_A.dg)**3),0)    
    
    print('----------------------------------------')
    print('                 Viga {:5.0f}           '.format(j+1))
    print('----------------------------------------')
    for i in np.arange(inicio,carga,passo):
        
        Met_A.esf(i,'distribuida')
        f_RRS,f_FMPS,f_EMAF,f_FMAV,f_FLT,ver_ELS = Met_A.verissimo()
        
        if  f_RRS  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por RRS')
            print('Posição de falha: {:5.2f}'.format(Met_A.pos_RRS))
            print('Flecha máxima = {:5.2f}'.format(Met_A.flecha))
            print('Flecha limite = {:5.2f}'.format(Met_A.flechalim))
            A_data[j,:] = np.array([i,'RRS',Met_A.pos_RRS,Met_A.flecha,      \
                                  Met_A.flechalim,Met_A.R_LeDG,Met_A.R_PeD0,\
                                      Met_A.R_DGeD0,Met_A.R_D0eD])
            break
        
        if  f_FMPS  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FMPS')
            print('Posição de falha: {:5.2f}'.format(Met_A.pos_FMPS))
            print('Flecha máxima = {:5.2f}'.format(Met_A.flecha))
            print('Flecha limite = {:5.2f}'.format(Met_A.flechalim))
            A_data[j,:] = np.array([i,'FMPS',Met_A.pos_FMPS,Met_A.flecha,    \
                                  Met_A.flechalim,Met_A.R_LeDG,Met_A.R_PeD0,\
                                      Met_A.R_DGeD0,Met_A.R_D0eD])
            break
    
        elif  f_EMAF  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por EMAF')
            print('Posição de falha: {:5.2f}'.format(Met_A.pos_EMAF))
            print('Flecha máxima = {:5.2f}'.format(Met_A.flecha))
            print('Flecha limite = {:5.2f}'.format(Met_A.flechalim))
            A_data[j,:] = np.array([i,'EMAF',Met_A.pos_EMAF,Met_A.flecha,      \
                                  Met_A.flechalim,Met_A.R_LeDG,Met_A.R_PeD0,\
                                      Met_A.R_DGeD0,Met_A.R_D0eD])
            break
        
        if  f_FMAV  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FMAV')
            print('Posição de falha: {:5.2f}'.format(Met_A.pos_FMAV))
            print('Flecha máxima = {:5.2f}'.format(Met_A.flecha))
            print('Flecha limite = {:5.2f}'.format(Met_A.flechalim))
            A_data[j,:] = np.array([i,'FMAV',Met_A.pos_FMAV,Met_A.flecha,    \
                                  Met_A.flechalim,Met_A.R_LeDG,Met_A.R_PeD0,\
                                      Met_A.R_DGeD0,Met_A.R_D0eD])
            break
        
        elif  f_FLT  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FLT')
            print('Posição de falha: {:5.2f}'.format(Met_A.pos_FLT))
            print('Flecha máxima = {:5.2f}'.format(Met_A.flecha))
            print('Flecha limite = {:5.2f}'.format(Met_A.flechalim))
            A_data[j,:] = np.array([i,'FLT',Met_A.pos_FLT,Met_A.flecha,      \
                                  Met_A.flechalim,Met_A.R_LeDG,Met_A.R_PeD0,\
                                      Met_A.R_DGeD0,Met_A.R_D0eD])
            break
    
        elif  ver_ELS  == 'flecha excedida':
            print('Carga de falha = {:5.2f}'.format(i)) 
            print('Falha por ELS')
            print('Flecha máxima = {:5.2f}'.format(Met_A.flecha))
            print('Flecha limite = {:5.2f}'.format(Met_A.flechalim))
            A_data[j,:] = np.array([i,'ELS','Centro',Met_A.flecha,          \
                                  Met_A.flechalim,Met_A.R_LeDG,Met_A.R_PeD0,\
                                      Met_A.R_DGeD0,Met_A.R_D0eD])
            break
    
    data_ver = pd.DataFrame(index = np.arange(n_vigas)+1,
                             data = A_data,
                             columns = ['Carga de falha','Tipo de falha',   \
                                        'Posição de falha','Flecha máx',    \
                                        'Flecha Limite','L/dg','p/d0','dg/d0',\
                                            'd0/d'])

print('----------------------------------------','\n')    

# Método AISC 31 (2016)

B_data = np.empty((n_vigas,9), dtype = 'object_')

print('----------------------------------------')
print('          Método ASIC (2016)            ')
print('----------------------------------------','\n')

for j in range(n_vigas):
    
    Met_B = viga_al(j)
    
    carga        = round(((30.72*2) * Met_A.melast * Met_A.Ix_O / (5*Met_A.dg)**3),1) 
    passo        = 0.1
    inicio       = round(((30.72/3) * Met_A.melast * Met_A.Ix_O / (20*Met_A.dg)**3),0)
    
    print('----------------------------------------')
    print('                 Viga {:5.0f}           '.format(j+1))
    print('----------------------------------------')
       
    for i in np.arange(inicio,carga,passo):
        
        Met_B.esf(i,'distribuida')
        f_RRS,f_FMVa,f_FMVm,f_FMV,f_FMA,f_FMAC,f_FMACnet,f_FLT,             \
        AISC_ELS=Met_B.AISC()
        
        if  f_RRS  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por RRS')
            print('Posição de falha: {:5.2f}'.format(Met_B.B_pos_RRS))
            print('Flecha máxima = {:5.2f}'.format(Met_B.fle_AISC))
            print('Flecha limite = {:5.2f}'.format(Met_B.flelim_AISC))
            B_data[j,:] = np.array([i,'RRS',Met_B.B_pos_RRS,Met_B.fle_AISC,    
                                  Met_B.flelim_AISC,Met_B.R_LeDG,Met_B.R_PeD0,\
                                      Met_B.R_DGeD0,Met_B.R_D0eD])
            break

        elif  f_FMVa  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FMV - Axial')
            print('Posição de falha: {:5.2f}'.format(Met_B.B_pos_FMVa))
            print('Flecha máxima = {:5.2f}'.format(Met_B.fle_AISC))
            print('Flecha limite = {:5.2f}'.format(Met_B.flelim_AISC))
            B_data[j,:] = np.array([i,'FMV_Axial',Met_B.B_pos_FMVa,           
                                  Met_B.fle_AISC,Met_B.flelim_AISC,\
                                      Met_B.R_LeDG,Met_B.R_PeD0,\
                                      Met_B.R_DGeD0,Met_B.R_D0eD])
            break

        elif  f_FMVm  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FMV - Flexão')
            print('Posição de falha: {:5.2f}'.format(Met_B.B_pos_FMVm))
            print('Flecha máxima = {:5.2f}'.format(Met_B.fle_AISC))
            print('Flecha limite = {:5.2f}'.format(Met_B.flelim_AISC))
            B_data[j,:] = np.array([i,'FMV_Flexão',Met_B.B_pos_FMVm,
                                  Met_B.fle_AISC,Met_B.flelim_AISC,
                                  Met_B.flelim_AISC,\
                                      Met_B.R_LeDG,Met_B.R_PeD0,\
                                      Met_B.R_DGeD0,Met_B.R_D0eD])
            break

        elif  f_FMV  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FMV')
            print('Posição de falha: {:5.2f}'.format(Met_B.B_pos_FMV))
            print('Flecha máxima = {:5.2f}'.format(Met_B.fle_AISC))
            print('Flecha limite = {:5.2f}'.format(Met_B.flelim_AISC))
            B_data[j,:] = np.array([i,'FMV',Met_B.B_pos_FMV,Met_B.fle_AISC,
                                  Met_B.flelim_AISC,Met_B.R_LeDG,Met_B.R_PeD0,\
                                      Met_B.R_DGeD0,Met_B.R_D0eD])
            break
        
        elif  f_FMA  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FMA') 
            print('Posição de falha: {:5.2f}'.format(Met_B.B_pos_FMA))
            print('Flecha máxima = {:5.2f}'.format(Met_B.fle_AISC))
            print('Flecha limite = {:5.2f}'.format(Met_B.flelim_AISC))
            B_data[j,:] = np.array([i,'FMA',Met_B.B_pos_FMA,Met_B.fle_AISC,
                                  Met_B.flelim_AISC,Met_B.R_LeDG,Met_B.R_PeD0,\
                                      Met_B.R_DGeD0,Met_B.R_D0eD])
            break

        elif  f_FMAC  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FMAC')
            print('Posição de falha: {:5.2f}'.format(Met_B.B_pos_FMC))
            print('Flecha máxima = {:5.2f}'.format(Met_B.fle_AISC))
            print('Flecha limite = {:5.2f}'.format(Met_B.flelim_AISC))
            B_data[j,:] = np.array([i,'FMA',Met_B.B_pos_FMC,Met_B.fle_AISC,
                                  Met_B.flelim_AISC,Met_B.R_LeDG,Met_B.R_PeD0,\
                                      Met_B.R_DGeD0,Met_B.R_D0eD])
            break

        elif  f_FMACnet  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por corte vertical na abertura')
            print('Posição de falha: {:5.2f}'.format(Met_B.B_pos_FMACnet))
            print('Flecha máxima = {:5.2f}'.format(Met_B.fle_AISC))
            print('Flecha limite = {:5.2f}'.format(Met_B.flelim_AISC))
            B_data[j,:] = np.array([i,'Falha por corte vertical na abertura',
                                  Met_B.B_pos_FMACnet,Met_B.fle_AISC,
                                  Met_B.flelim_AISC,Met_B.R_LeDG,Met_B.R_PeD0,\
                                      Met_B.R_DGeD0,Met_B.R_D0eD])
            break

        elif  f_FLT  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FLT')
            print('Início de falha: {:5.2f}'.format(Met_B.B_pos_FLT))
            print('Flecha máxima = {:5.2f}'.format(Met_B.fle_AISC))
            print('Flecha limite = {:5.2f}'.format(Met_B.flelim_AISC))
            B_data[j,:] = np.array([i,'FLT',Met_B.B_pos_FLT,
                                  Met_B.fle_AISC,Met_B.flelim_AISC,\
                                      Met_B.flelim_AISC,\
                                      Met_B.R_LeDG,Met_B.R_PeD0,\
                                      Met_B.R_DGeD0,Met_B.R_D0eD])
            break
        
        elif  AISC_ELS  == 'flecha excedida':
            print('Carga de falha = {:5.2f}'.format(i)) 
            print('Falha por ELS')
            print('Flecha máxima = {:5.2f}'.format(Met_B.fle_AISC))
            print('Flecha limite = {:5.2f}'.format(Met_B.flelim_AISC))
            B_data[j,:] = np.array([i,'ELS','Centro',Met_B.fle_AISC,
                                  Met_B.flelim_AISC,\
                                      Met_B.R_LeDG,Met_B.R_PeD0,\
                                      Met_B.R_DGeD0,Met_B.R_D0eD])
            break
        
        elif  Met_B.Lim_AISC  == 'false': 
            print('Viga fora do limite')
            B_data[j,:] = np.array(['Viga fora do limite','Viga fora do limite',
                                  'Viga fora do limite','Viga fora do limite',
                                  'Viga fora do limite','Viga fora do limite',\
                                  'Viga fora do limite','Viga fora do limite',\
                                  'Viga fora do limite'])
            break  
    
    data_AISC = pd.DataFrame(index = np.arange(n_vigas)+1,
                             data = B_data,
                             columns = ['Carga de falha','Tipo de falha',   \
                                        'Posição de falha','Flecha máx',    \
                                        'Flecha Limite','L/dg','p/d0','dg/d0',\
                                            'd0/d'])    
        
    print('Classificação de seção:'.format(Met_B.comp))
    
print('----------------------------------------','\n')    

# Método SCI No 100 and BS5950 

C_data = np.empty((n_vigas,9), dtype = 'object_')

print('----------------------------------------')
print('     Método SCI No 100 and BS5950       ')
print('----------------------------------------','\n')

for j in range(n_vigas):
    
    Met_C = viga_al(j)
    
    carga        = round(((30.72*2) * Met_A.melast * Met_A.Ix_O / (5*Met_A.dg)**3),1) 
    passo        = 0.1
    inicio       = round(((30.72/3) * Met_A.melast * Met_A.Ix_O / (20*Met_A.dg)**3),0)
    
    print('----------------------------------------')
    print('                 Viga {:5.0f}           '.format(j+1))
    print('----------------------------------------')
        
    for i in np.arange(inicio,carga,passo):
        
        Met_C.esf(i,'distribuida')
        f_FMPS,f_RRS,f_VCV,f_FMA,f_FMV,f_FLT,ver_ELS=Met_C.SCI()
        
        if  f_FMPS  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FMPS')
            print('Posição de falha: {:5.2f}'.format(Met_C.C_pos_FMPS))
            print('Flecha máxima = {:5.2f}'.format(Met_C.SCI_Umax))
            print('Flecha limite = {:5.2f}'.format(Met_C.SCI_Ulim))
            C_data[j,:] = np.array([i,'FMPS',Met_C.C_pos_FMPS,Met_C.SCI_Umax,
                                 Met_C.SCI_Ulim,Met_C.R_LeDG,Met_C.R_PeD0,\
                                     Met_C.R_DGeD0,Met_C.R_D0eD])
            break

        elif  f_RRS  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por RRS')
            print('Posição de falha: {:5.2f}'.format(Met_C.C_pos_RRS))
            print('Flecha máxima = {:5.2f}'.format(Met_C.SCI_Umax))
            print('Flecha limite = {:5.2f}'.format(Met_C.SCI_Ulim))
            C_data[j,:] = np.array([i,'RRS',Met_C.C_pos_RRS,Met_C.SCI_Umax,
                                 Met_C.SCI_Ulim,Met_C.R_LeDG,Met_C.R_PeD0,\
                                     Met_C.R_DGeD0,Met_C.R_D0eD])
            break
        
        elif  f_VCV  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por VCV')
            print('Posição de falha: {:5.2f}'.format(Met_C.C_pos_VCV))
            print('Flecha máxima = {:5.2f}'.format(Met_C.SCI_Umax))
            print('Flecha limite = {:5.2f}'.format(Met_C.SCI_Ulim))
            C_data[j,:] = np.array([i,'VCV',Met_C.C_pos_VCV,Met_C.SCI_Umax,
                                 Met_C.SCI_Ulim,Met_C.R_LeDG,Met_C.R_PeD0,\
                                     Met_C.R_DGeD0,Met_C.R_D0eD])
            break
        
        elif  f_FMA  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FMA')
            print('Posição de falha: {:5.2f}'.format(Met_C.C_pos_FMA))
            print('Flecha máxima = {:5.2f}'.format(Met_C.SCI_Umax))
            print('Flecha limite = {:5.2f}'.format(Met_C.SCI_Ulim))
            C_data[j,:] = np.array([i,'FMA',Met_C.C_pos_FMA,Met_C.SCI_Umax,
                                  Met_C.SCI_Ulim,Met_C.R_LeDG,Met_C.R_PeD0,\
                                      Met_C.R_DGeD0,Met_C.R_D0eD])
            break
        
        elif  f_FMV  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FMV')
            print('Posição de falha: {:5.2f}'.format(Met_C.C_pos_FMV))
            print('Flecha máxima = {:5.2f}'.format(Met_C.SCI_Umax))
            print('Flecha limite = {:5.2f}'.format(Met_C.SCI_Ulim))
            C_data[j,:] = np.array([i,'FMV',Met_C.C_pos_FMV,Met_C.SCI_Umax,
                                 Met_C.SCI_Ulim,Met_C.R_LeDG,Met_C.R_PeD0,\
                                     Met_C.R_DGeD0,Met_C.R_D0eD])
            break
        
        elif  f_FLT  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FLT')
            print('Posição de falha: {:5.2f}'.format(Met_C.C_pos_FLT))
            print('Flecha máxima = {:5.2f}'.format(Met_C.SCI_Umax))
            print('Flecha limite = {:5.2f}'.format(Met_C.SCI_Ulim))
            C_data[j,:] = np.array([i,'FLT',Met_C.C_pos_FLT,Met_C.SCI_Umax,
                                 Met_C.SCI_Ulim,Met_C.R_LeDG,Met_C.R_PeD0,\
                                     Met_C.R_DGeD0,Met_C.R_D0eD])
            break

        elif  ver_ELS  == 'flecha excedida':
            print('Carga de falha = {:5.2f}'.format(i)) 
            print('Falha por ELS - Verificação simplificada')
            print('Flecha máxima = {:5.2f}'.format(Met_C.SCI_Umax))
            print('Flecha limite = {:5.2f}'.format(Met_C.SCI_Ulim))
            C_data[j,:] = np.array([i,'ELS','Centro',Met_C.SCI_Umax,
                                  Met_C.SCI_Ulim,Met_C.R_LeDG,Met_C.R_PeD0,\
                                      Met_C.R_DGeD0,Met_C.R_D0eD])
            break
        
        elif  Met_C.Lim_SCI  == 'false': 
            print('Viga fora do limite')
            C_data[j,:] = np.array(['Viga fora do limite','Viga fora do limite',
                                  'Viga fora do limite','Viga fora do limite',
                                  'Viga fora do limite','Viga fora do limite',\
                                  'Viga fora do limite','Viga fora do limite',\
                                  'Viga fora do limite'])
            break
    
    data_SCI = pd.DataFrame(index = np.arange(n_vigas)+1,
                              data = C_data,
                              columns = ['Carga de falha','Tipo de falha',   \
                                         'Posição de falha','Flecha máx',    \
                                         'Flecha Limite','L/dg','p/d0','dg/d0',\
                                             'd0/d'])
print('----------------------------------------','\n')
    
# Método GRILO et al. (2018)

D_data = np.empty((n_vigas,9), dtype = 'object_')

print('----------------------------------------')
print('       Método GRILO et al. (2018)       ')
print('----------------------------------------','\n')

for j in range(n_vigas):
    
    Met_D = viga_al(j)
    
    carga        = round(((30.72*2) * Met_A.melast * Met_A.Ix_O / (5*Met_A.dg)**3),1) 
    passo        = 0.1
    inicio       = round(((30.72/3) * Met_A.melast * Met_A.Ix_O / (20*Met_A.dg)**3),0)
    
    print('----------------------------------------')
    print('                 Viga {:5.0f}           '.format(j+1))
    print('----------------------------------------')
        
    for i in np.arange(inicio,carga,passo):
        Met_D.esf(i,'distribuida')
        f_lim,f_FMA,D_ver_ELS,f_FMPS,f_RRS,f_FLT=Met_D.grilo()
        
        if  f_FMA  == 'Falha por FMA':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FMAV')
            print('Posição de falha: {:5.2f}'.format(Met_D.D_pos_FMA))
            print('Flecha máxima = {:5.2f}'.format(Met_D.fle_grilo))
            print('Flecha limite = {:5.2f}'.format(Met_D.flelim_grilo))
            D_data[j,:] = np.array([i,'FMA',Met_D.D_pos_FMA,Met_D.fle_grilo,
                                 Met_D.flelim_grilo,Met_D.R_LeDG,Met_D.R_PeD0,\
                                     Met_D.R_DGeD0,Met_D.R_D0eD])
            break

        elif  D_ver_ELS  == 'flecha excedida':
            print('Carga de falha = {:5.2f}'.format(i)) 
            print('Falha por ELS')
            print('Flecha máxima = {:5.2f}'.format(Met_D.fle_grilo))
            print('Flecha limite = {:5.2f}'.format(Met_D.flelim_grilo))
            D_data[j,:] = np.array([i,'ELS','Centro',Met_D.fle_grilo,
                                  Met_D.flelim_grilo,Met_D.R_LeDG,Met_D.R_PeD0,\
                                      Met_D.R_DGeD0,Met_D.R_D0eD])
            break
            
        elif  f_lim  == 'false': 
            print('Viga fora do limite')
            D_data[j,:] = np.array(['Viga fora do limite','Viga fora do limite',
                                  'Viga fora do limite','Viga fora do limite',
                                  'Viga fora do limite','Viga fora do limite',\
                                  'Viga fora do limite','Viga fora do limite',\
                                  'Viga fora do limite'])
            break
 
# Copiado do método Veríssimo et al. (2012)   

        elif  f_RRS  == 'falha':
             print('Carga de falha = {:5.2f}'.format(i))
             print('Falha por RRS')
             print('Posição de falha: {:5.2f}'.format(Met_D.pos_RRS))
             print('Flecha máxima = {:5.2f}'.format(Met_D.fle_grilo))
             print('Flecha limite = {:5.2f}'.format(Met_D.flelim_grilo))
             D_data[j,:] = np.array([i,'RRS',Met_D.pos_RRS,Met_D.fle_grilo,   \
                                  Met_D.flelim_grilo,Met_D.R_LeDG,Met_D.R_PeD0,\
                                      Met_D.R_DGeD0,Met_D.R_D0eD])
             break
         
        elif  f_FMPS  == 'falha':
             print('Carga de falha = {:5.2f}'.format(i))
             print('Falha por FMPS')
             print('Posição de falha: {:5.2f}'.format(Met_D.pos_FMPS))
             print('Flecha máxima = {:5.2f}'.format(Met_D.fle_grilo))
             print('Flecha limite = {:5.2f}'.format(Met_D.flelim_grilo))
             D_data[j,:] = np.array([i,'FMPS',Met_D.pos_FMPS,Met_D.fle_grilo, \
                                  Met_D.flelim_grilo,Met_D.R_LeDG,Met_D.R_PeD0,\
                                      Met_D.R_DGeD0,Met_D.R_D0eD])
             break
         
        elif  f_FLT  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FLT')
            print('Posição de falha: {:5.2f}'.format(Met_D.pos_FLT))
            print('Flecha máxima = {:5.2f}'.format(Met_D.fle_grilo))
            print('Flecha limite = {:5.2f}'.format(Met_D.flelim_grilo))
            D_data[j,:] = np.array([i,'FLT',Met_D.pos_FLT,Met_D.fle_grilo,    \
                                 Met_D.flelim_grilo,Met_D.R_LeDG,Met_D.R_PeD0,\
                                     Met_D.R_DGeD0,Met_D.R_D0eD])
            break               
     
    data_grilo = pd.DataFrame(index = np.arange(n_vigas)+1,
                              data = D_data,
                              columns = ['Carga de falha','Tipo de falha',   \
                                         'Posição de falha','Flecha máx',    \
                                         'Flecha Limite','L/dg','p/d0','dg/d0',\
                                             'd0/d'])

print('----------------------------------------','\n')

# Anexo N / Eurocode 3 (ENV 1993-1-1:1992/A2:1998) 

E_data = np.empty((n_vigas,9), dtype = 'object_')

print('----------------------------------------')
print('          Anexo N / Eurocode 3          ')
print('----------------------------------------','\n')

for j in range(n_vigas):
    
    Met_E = viga_al(j)
    
    carga        = round(((30.72*2) * Met_A.melast * Met_A.Ix_O / (5*Met_A.dg)**3),1) 
    passo        = 0.1
    inicio       = round(((30.72/3) * Met_A.melast * Met_A.Ix_O / (20*Met_A.dg)**3),0)
    
    print('----------------------------------------')
    print('                 Viga {:5.0f}           '.format(j+1))
    print('----------------------------------------')
        
    for i in np.arange(inicio,carga,passo):
                                              
        Met_E.esf(i,'distribuida')
        f_FMPS,f_FMV,f_FMA,f_VCV,f_RRS,ver_ELS=Met_E.annexN()
        
        if  f_FMPS  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FMPS')
            print('Posição de falha: {:5.2f}'.format(Met_E.E_pos_FMPS))
            print('Flecha máxima = {:5.2f}'.format(Met_E.annexN_Umax))
            print('Flecha limite = {:5.2f}'.format(Met_E.annexN_Ulim))
            E_data[j,:] = np.array([i,'FMPS',Met_E.E_pos_FMPS,Met_E.annexN_Umax,
                                 Met_E.annexN_Ulim,Met_E.R_LeDG,Met_E.R_PeD0,\
                                     Met_E.R_DGeD0,Met_E.R_D0eD])
            break
        
        elif  f_FMV  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FMV')
            print('Posição de falha: {:5.2f}'.format(Met_E.E_pos_FMV))
            print('Flecha máxima = {:5.2f}'.format(Met_E.annexN_Umax))
            print('Flecha limite = {:5.2f}'.format(Met_E.annexN_Ulim))
            E_data[j,:] = np.array([i,'FMV',Met_E.E_pos_FMV,Met_E.annexN_Umax,
                                 Met_E.annexN_Ulim,Met_E.R_LeDG,Met_E.R_PeD0,\
                                     Met_E.R_DGeD0,Met_E.R_D0eD])
            break
        
        elif  f_FMA  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FMA')
            print('Posição de falha: {:5.2f}'.format(Met_E.E_pos_FMA))
            print('Flecha máxima = {:5.2f}'.format(Met_E.annexN_Umax))
            print('Flecha limite = {:5.2f}'.format(Met_E.annexN_Ulim))
            E_data[j,:] = np.array([i,'FMA',Met_E.E_pos_FMA,Met_E.annexN_Umax,
                                 Met_E.annexN_Ulim,Met_E.R_LeDG,Met_E.R_PeD0,\
                                     Met_E.R_DGeD0,Met_E.R_D0eD])
            break

        elif  f_RRS  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por RRS')
            print('Posição de falha: {:5.2f}'.format(Met_E.E_pos_RRS))
            print('Flecha máxima = {:5.2f}'.format(Met_E.annexN_Umax))
            print('Flecha limite = {:5.2f}'.format(Met_E.annexN_Ulim))
            E_data[j,:] = np.array([i,'RRS',Met_E.E_pos_RRS,Met_E.annexN_Umax,
                                 Met_E.annexN_Ulim,Met_E.R_LeDG,Met_E.R_PeD0,\
                                     Met_E.R_DGeD0,Met_E.R_D0eD])
            break
        
        elif  f_VCV  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por VCV')
            print('Posição de falha: {:5.2f}'.format(Met_E.E_pos_VCV))
            print('Flecha máxima = {:5.2f}'.format(Met_E.annexN_Umax))
            print('Flecha limite = {:5.2f}'.format(Met_E.annexN_Ulim))
            E_data[j,:] = np.array([i,'VCV',Met_E.E_pos_VCV,Met_E.annexN_Umax,
                                 Met_E.annexN_Ulim,Met_E.R_LeDG,Met_E.R_PeD0,\
                                     Met_E.R_DGeD0,Met_E.R_D0eD])
            break
        
        # elif  f_FLT  == 'falha':
        #     print('Carga de falha = {:5.2f}'.format(i))
        #     print('Falha por FLT')
        #     print('Posição de falha: {:5.2f}'.format(Met_E.E_pos_FLT))
        #     print('Flecha máxima = {:5.2f}'.format(Met_E.annexN_Umax))
        #     print('Flecha limite = {:5.2f}'.format(Met_E.annexN_Ulim))
        #     pla[j,:] = np.array([i,'FLT',Met_E.E_pos_FMPS,Met_E.annexN_Umax,
        #                    Met_E.annexN_Ulim,Met_E.R_LeDG,Met_E.R_PeD0,\
        #                    Met_E.R_DGeD0,Met_E.R_D0eD])
        #     break

# Comentar qual verificação do ELS quer fazer para esse método
        
        elif  ver_ELS  == 'flecha excedida':
            print('Carga de falha = {:5.2f}'.format(i)) 
            print('Falha por ELS')
            print('Flecha máxima = {:5.2f}'.format(Met_E.annexN_Umax))
            print('Flecha limite = {:5.2f}'.format(Met_E.annexN_Ulim))
            E_data[j,:] = np.array([i,'ELS','Centro',Met_E.annexN_Umax,
                                  Met_E.annexN_Ulim,Met_E.R_LeDG,Met_E.R_PeD0,\
                                      Met_E.R_DGeD0,Met_E.R_D0eD])
            break
        
        elif  Met_E.Lim_annexN  == 'false': 
            print('Viga fora do limite')
            E_data[j,:] = np.array(['Viga fora do limite','Viga fora do limite',
                                  'Viga fora do limite','Viga fora do limite',
                                  'Viga fora do limite','Viga fora do limite',\
                                  'Viga fora do limite','Viga fora do limite',\
                                  'Viga fora do limite'])
            break
    
    data_annexN = pd.DataFrame(index = np.arange(n_vigas)+1,
                              data = E_data,
                              columns = ['Carga de falha','Tipo de falha',   \
                                         'Posição de falha','Flecha máx',    \
                                         'Flecha Limite','L/dg','p/d0','dg/d0',\
                                             'd0/d'])

print('----------------------------------------','\n')

print('Salvando em planilha Excel...')
writer = pd.ExcelWriter('Results_AlveolarBeams.xlsx')
data_ver.to_excel(writer,'Veríssimo et al. (2012)',index = 'false')
data_AISC.to_excel(writer,'AISC 31 (2016)',index = 'false')
data_SCI.to_excel(writer,'SCI No 100 and BS5950',index = 'false')
data_grilo.to_excel(writer,'GRILO et al. (2018)',index = 'false')
data_annexN.to_excel(writer,'Anexo N (Eurocode 3)',index = 'false')
   
# writer.save()   # Descontinuado
writer.close()

print()

# Verificação das vigas originais

O_data = np.empty((n_vigas,5), dtype = 'object_')

print('----------------------------------------')
print('       Perfis originais (NBR 8800)      ')
print('----------------------------------------','\n')

for j in range(n_vigas):
    
    teste = viga_al(j)
    
    carga        = round(((30.72*2) * Met_A.melast * Met_A.Ix_O / (5*Met_A.d)**3),1) 
    passo        = 0.1
    inicio       = round(((30.72/3) * Met_A.melast * Met_A.Ix_O / (20*Met_A.d)**3),0)
    
    print('----------------------------------------')
    print('                 Viga {:5.0f}           '.format(j+1))
    print('----------------------------------------')
    for i in np.arange(inicio,carga,passo):
        
        teste.esf(i,'distribuida')
        f_FLM,f_FLA,f_FCR,f_FLT,f_ELS = teste.original()
        
        if  f_FLM  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FLM')
            print('Posição de falha: {:5.2f}'.format(teste.O_pos_FLM))
            print('Flecha máxima = {:5.2f}'.format(teste.Ori_fle))
            print('Flecha limite = {:5.2f}'.format(teste.Ori_flelim))
            O_data[j,:] = np.array([i,'FLM',teste.O_pos_FLM,teste.Ori_fle,
                                  teste.Ori_flelim])
            break
        
        if  f_FLA  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FLA')
            print('Posição de falha: {:5.2f}'.format(teste.O_pos_FLA))
            print('Flecha máxima = {:5.2f}'.format(teste.Ori_fle))
            print('Flecha limite = {:5.2f}'.format(teste.Ori_flelim))
            O_data[j,:] = np.array([i,'FLA',teste.O_pos_FLA,teste.Ori_fle,
                                 teste.Ori_flelim])
            break

        elif  f_FCR  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FCR')
            print('Posição de falha: {:5.2f}'.format(teste.O_pos_FCR))
            print('Flecha máxima = {:5.2f}'.format(teste.Ori_fle))
            print('Flecha limite = {:5.2f}'.format(teste.Ori_flelim))
            O_data[j,:] = np.array([i,'FCR',teste.O_pos_FCR,teste.Ori_fle,
                                 teste.Ori_flelim])
            break
        
        elif  f_FLT  == 'falha':
            print('Carga de falha = {:5.2f}'.format(i))
            print('Falha por FLT')
            print('Posição de falha: {:5.2f}'.format(teste.O_pos_FLT))
            print('Flecha máxima = {:5.2f}'.format(teste.Ori_fle))
            print('Flecha limite = {:5.2f}'.format(teste.Ori_flelim))
            O_data[j,:] = np.array([i,'FLT',teste.O_pos_FLT,teste.Ori_fle,
                                 teste.Ori_flelim])
            break
        
        elif  f_ELS  == 'flecha excedida':
            print('Carga de falha = {:5.2f}'.format(i)) 
            print('Falha por ELS')
            print('Flecha máxima = {:5.2f}'.format(teste.Ori_fle))
            print('Flecha limite = {:5.2f}'.format(teste.Ori_flelim))
            O_data[j,:] = np.array([i,'ELS','Centro',teste.Ori_fle,
                                 teste.Ori_flelim])
            break

    data_ori = pd.DataFrame(index = np.arange(n_vigas)+1,
                              data = O_data,
                              columns = ['Carga de falha','Tipo de falha',
                                         'Posição de falha','Flecha máx',
                                         'Flecha Limite'])

print('----------------------------------------','\n')
 

print('Salvando em planilha Excel...')
writer = pd.ExcelWriter('Results_OriginalISection.xlsx')
data_ori.to_excel(writer,'NBR 8800',index = 'false')

# writer.save()
writer.close() 