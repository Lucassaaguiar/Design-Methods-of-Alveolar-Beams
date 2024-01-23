##############################################################################
#               Programa de cálculo de vigas alveolares                      #
#                                                                            #
# Objetivo: Realizar leitura dos dados de entrada, verificações e            # 
#           dimensionamento dos métodos para vigas alveolares                #
# Programaddor: Lucas Alves de Aguiar                                        #
# Início: 15/09/2022                                                         #
# Revisão: 29/01/2023                                                        #  
# Atualização: 22/01/2024                                                    # 
##############################################################################

import numpy as np
import pandas as pd
from scipy.interpolate import interp2d

#----------------------------------------------------------------------------#
# 1 - Inicializando a Classe para vigas celulares                            #
#----------------------------------------------------------------------------#

class viga_al():
    '''
    Classe das vigas alveolares, contém as propriedades de dados de 
    entrda do sistema das vigas do tipo celulares.
    '''
    
    def __init__(self,nl):
        '''
        Criando o objeto: Vigas - contendo os dados de entrada (de)
        '''
        de    = pd.read_excel('exemples.xlsx', sheet_name='dados')              # Importando dados de planilha de Excel
        de    = de.values
        dm    = pd.read_excel('exemples.xlsx', sheet_name='materiais')          # Importando dados de planilha de Excel
        dm    = dm.values
        
        # Parametros para grilo
        self.op_a    = pd.read_excel('coefficients_griloFMA.xlsx', 
                       sheet_name='a',index_col =0)
        self.op_b    = pd.read_excel('coefficients_griloFMA.xlsx', 
                       sheet_name='b',index_col =0)
        self.op_c    = pd.read_excel('coefficients_griloFMA.xlsx', 
                       sheet_name='c',index_col =0)
        self.op_d    = pd.read_excel('coefficients_griloFMA.xlsx', 
                       sheet_name='d',index_col =0)
        self.op_e    = pd.read_excel('coefficients_griloFMA.xlsx', 
                       sheet_name='e',index_col =0)
        
        # Leitura de dados
        
        self.dg     = de[nl,0]    
        self.d      = de[nl,1]
        self.d0     = de[nl,2]
        self.h0     = de[nl,3]
        self.l      = de[nl,4]
        self.bwe    = de[nl,5]
        self.bwmin  = de[nl,6]
        self.bw     = de[nl,7]
        self.p      = de[nl,8]
        self.n      = int(de[nl,9])
        self.lb     = de[nl,10]
        self.hp     = de[nl,11]
        self.bf     = de[nl,12]
        self.tf     = de[nl,13]
        self.tw     = de[nl,14]
        self.bi     = de[nl,15]
        self.coefseg= de[nl,16]
        self.sist   = de[nl,17]
                       
#----------------------------------------------------------------------------#
# 1.1 - Cálculo das Propriedades Geométricas                                 #
#----------------------------------------------------------------------------#

# Propriedades da seção transversal do T
 
        ht      = (self.dg-self.d0)/2                                           # Altura do T                                                                                     
        At      = self.bf*self.tf + (ht-self.tf)*self.tw                        # Área do T
        Anet    = 2*At                                                          # Área dos T's
        yl      = ((self.bf*self.tf**2) + (self.tw*ht**2) - 
                  (self.tw*self.tf**2)) / (2*(self.bf*self.tf + ht*self.tw - 
                   self.tf*self.tw))                                            # centro geométrico do T em y (Eixo de corte x)
        Ix_t    = (((self.bf*self.tf**3)/12) + (self.bf*self.tf*
                  (yl-self.tf/2)**2) + ((self.tw*(ht-self.tf)**3)/12) +      
                  (self.tw*(ht-self.tf)*(yl-(ht+self.tf)/2)**2))                # Inércia x do T
        y0      = self.d0/2 + ht - yl                                           # Distância da linha neutra do T ao eixo da viga
        ya      = ht - yl                                                       # Distância da linha neutra do T ao todo do alveolo
        yb      = yl - self.tf                                                  # Distância da linha neutra do T até o início da flange   
        deffec  = self.dg-2*yl                                                  # Distância entre as linhas neutras dos T
        Iy_t    = (self.tf*self.bf**3)/12 + ((ht-self.tf)*self.tw**3)/12        # Inércia (y) da seção T
        rx_t    = np.sqrt(Ix_t/At)                                              # Raio de giração em relação ao eixo x 
        ry_t    = np.sqrt(Iy_t/At)                                              # Raio de giração em relação ao eixo y
        Cw_t    = (1/36) * (((self.bf**3*self.tf**3)/4) + 
                  ((ht-self.tf/2)**3) * (self.tw**3))                           # Constante de empenhamento da seção transversal do T  
        J_t     = (2/3 * self.bf*self.tf**3 + (ht-self.tf)*self.tw**3/3)        # Constante de rotação da seção transversal do T
        
        if At/(2*self.bf) < self.tf:
            Zx_t= (((self.tw*ht**2)/2) + ((self.bf*self.tf**2)/4) - \
                   ((ht*self.tf *self.tw)/2) - (((ht-self.tf)**2)/(4*self.bf)))
        else:
            Zx_t= ((self.tw*(ht-self.tf)**2/4) + (self.bf*ht*self.tf/2) - \
                  ((self.bf**2*self.tf**2)/(4*self.tw)))                        # Módulo da resistência plástica X (Horizontal)
     
        Wx_tc   = Ix_t/ya                                                       # Módulo de resistência elástica x da seção T da região comprimida  
        Wx_tt   = Iy_t/yl                                                       # Módulo de resistência elástica x da seção T da região tracionada
        Wy_t    = 2*Iy_t/self.bf                                                # Módulo de resistência elástica y da seção T                 
        if self.bi != 0 :                                                       # Verificação de viga castelada (se tiver bi = castelada)
           Hexp = (self.h0+self.hp)/2                                           # Altura de metáde da abertura 
           a0   = 2*self.bi + self.bw                                           # Largura do Alveolo  
        else:
            
           Hexp = self.d0/2                                                     # Altura de metáde da abertura
           a0   = self.d0                                                       # Largura do Alveolo
        
# Propriedades da seção transversal sem alveolo (alma cheia)

        Ag      = self.bf*self.tf*2 + self.tw*(self.dg-2*self.tf)               # Área da seção do montante da alma cheia        
        Ix      = 2*(self.bf*self.tf**3/12 + self.bf*self.tf*((self.dg/2)-
                  (self.tf/2))**2) + ((self.dg-2*self.tf)**3*self.tw)/12        # Inércia (x) da seção de alma cheia
        Wx      = 2*Ix/self.dg                                                  # Módulo de resistência elástica da seção de alma cheia (x)
        rx      = np.sqrt(Ix/Ag)                                                # Raio de giração em relação ao eixo x
        Zx      = (self.bf*self.tf*(self.dg-self.tf) + 
                  (((self.dg-2*self.tf)**2)*self.tw)/4)                         # Módulo de resistência plástica x
        Iy      = (2*(self.tf*self.bf**3)/12 + 
                  ((self.dg-2*self.tf)*self.tw**3/12))                          # Inércia (y) da seção de alma cheia
        Wy      = 2*Iy/self.bf                                                  # Módulo de resistência elástica da seção de alma cheia (y)
        ry      = np.sqrt(Iy/Ag)                                                # Raio de giração em relação ao eixo y (Horizontal)
        Cw      = Iy*(self.dg-self.tf)**2/4                                     # Constante de empenhamento da seção transversal
        J       = 2/3 * self.bf*self.tf**3 + (self.dg-self.tf)*self.tw**3/3     # Constante de rotação da seção transversal com alma cheia  
           
# Propriedades da seção transversal com alveolo (alma vazada) - 
# Faria et al. (2020)

        Ix_vz   = (self.bf*self.tf**3)/6  + (self.bf*self.tf*
                  (self.dg-self.tf)**2)/2 + ((self.dg-2*self.tf)**3-
                   self.d0**3)*self.tw/12                                       # Inércia (x) da seção vazada da viga
        Wx_vz   = 2*Ix_vz/self.dg                                               # Módulo da resistência elástica da seção (Eixo X) / Sxv
        rx_vz   = np.sqrt(Ix_vz/Anet)                                           # Raio de giração em relação ao eixo X (Horizontal)
        Zx_vz   = self.bf*self.tf*(self.dg-self.tf) + ((self.dg-2*self.tf)**2-
                  self.d0**2)*self.tw/4                                         # Módulo da resistência plástica X (Horizontal) 
        Iy_vz   = self.tf*self.bf**3/6 + (self.dg-2*self.tf-
                  self.d0)*self.tw**3/12                                        # Inércia (y) da seção vazada da viga
        Wy_vz   = 2*Iy_vz/self.bf                                               # Módulo da resistência elástica da seção (Eixo X) / Sxv
        ry_vz   = np.sqrt(Iy_vz/Anet)                                           # Raio de giração em relação ao eixo Y (Vertica
        Cw_vz   = (self.bf**3)*((self.dg-self.tf)**2)*self.tf/24                # Constante de empenhamento da seção Transversal com alma vazada
        J_vz    = J - (self.d0/3)*self.tw**3                                    # Constante de rotação da seção transversal com alma vazada
        Ag_vz   = Ag - (self.d0*self.tw)                                        # Área da seção vazada da viga (região do alveolo)

# Propriedades da seção transversal da viga original

        self.Ag_O    = self.bf*self.tf*2 + self.tw*(self.d-2*self.tf)         
        self.Ix_O    = 2*(self.bf*self.tf**3/12 + self.bf*self.tf*((self.d/2)-
                       (self.tf/2))**2) + ((self.d-2*self.tf)**3*self.tw)/12    # Inércia (x) da seção do perfil original da viga
        self.Wx_O    = 2*self.Ix_O/self.d                                       # Módulo de resistência elástica da seção do perfil original (x)
        self.rx_O    = np.sqrt(self.Ix_O/self.Ag_O)                             # Raio de giração em relação ao eixo x
        self.Zx_O    = (self.bf*self.tf*(self.d-self.tf) + 
                       ((self.d-2*self.tf)**2)*self.tw/4)                       # Módulo de resistência plástica x
        self.Iy_O    = (2*(self.tf*self.bf**3)/12 + 
                       ((self.d-2*self.tf)*self.tw**3/12))                      # Inércia (y) da seção do perfil original
        self.Wy_O    = 2*self.Iy_O/self.bf                                      # Módulo de resistência elástica da seção original (y)
        self.ry_O    = np.sqrt(self.Iy_O/self.Ag_O)                             # Raio de giração em relação ao eixo y (Horizontal)
        # Zy_O   = 
        self.Cw_O    = self.Iy_O*(self.d-self.tf)**2/4                          # Constante de empenhamento da seção Transversal do perfil original
        self.J_O     = (2/3 * self.bf*self.tf**3 + 
                       (self.d-self.tf)*self.tw**3)/3                           # Constante de rotação da seção transversal do perfil original  
        self.Ag_O    = self.bf*self.tf*2 + self.tw*(self.d-2*self.tf)           # Área da seção do perfil original

# Seção T crítica 

        y_crit       = np.sqrt((0.5*self.d0)**2 - (0.225*self.d0)**2)
        ht_crit      = self.d0/2 - y_crit + ht
        yl_crit      = ((self.bf*self.tf**2) + (self.tw*ht_crit**2) - 
                        (self.tw*self.tf**2)) / (2*(self.bf*self.tf + 
                         ht_crit*self.tw - self.tf*self.tw))  
        deffec_crit  = self.dg - 2 * yl_crit
        At_crit      = (ht_crit-self.tf)*self.tw + self.bf * self.tf       
        Ix_t_crit    = (((self.bf*self.tf**3)/12)+(self.bf*self.tf*
                        (yl_crit-self.tf/2)**2)+((self.tw*(ht_crit-self.tf)**3)
                         /12)+(self.tw*(ht_crit-self.tf)*
                         (yl_crit-(ht_crit+self.tf)/2)**2))
        Iy_t_crit    = (self.tf*self.bf**3)/12 + ((ht_crit-
                        self.tf)*self.tw**3)/12           
        rx_t_crit    = np.sqrt(Ix_t_crit/At_crit) 
        ry_t_crit    = np.sqrt(Iy_t_crit/At_crit)       
        Cw_t_crit    = (1/36) * (((self.bf**3*self.tf**3)/4) + 
                       ((ht_crit-self.tf/2)**3) * (self.tw**3))
        
        if At_crit/(2*self.bf) < self.tf:
            Zx_t_crit= (((self.tw*ht_crit**2)/2) + ((self.bf*self.tf**2)/4) - \
                   ((ht_crit*self.tf *self.tw)/2) - (((ht_crit-self.tf)**2)/(4*self.bf)))
        else:
           Zx_t_crit= ((self.tw*(ht_crit-self.tf)**2/4) + (self.bf*ht_crit*self.tf/2) - \
                  ((self.bf**2*self.tf**2)/(4*self.tw)))
        
        J_t_crit     = (2/3 * self.bf*self.tf**3 + (ht_crit-
                        self.tf)*self.tw**3/3)   
        Wx_tc_crit   = Ix_t_crit / (ht_crit - yl_crit)
        Wx_tt_crit   = Iy_t_crit / yl_crit
        Wy_t_crit    = 2*Iy_t_crit / self.bf                        
        Ix_vz_crit   = 2 * Ix_t_crit + 2 * At_crit*(deffec_crit/2)**2                           
        Wx_vz_crit   = Ix_vz_crit/(self.dg/2)
        Zx_vz_crit   = 2 * At_crit * deffec_crit / 2
        
#----------------------------------------------------------------------------#
# 1.2 - Associando as propriedades geométria à classe de vigas alveoladas    #
#----------------------------------------------------------------------------#   

        self.ht    = ht
        self.At    = At
        self.Anet  = Anet
        self.yl    = yl
        self.Ix_t  = Ix_t
        self.y0    = y0
        self.ya    = ya
        self.yb    = yb
        self.deffec= deffec
        self.Iy_t  = Iy_t
        self.rx_t  = rx_t
        self.ry_t  = ry_t
        self.Cw_t  = Cw_t
        self.J_t   = J_t
        self.Zx_t  = Zx_t
        self.Wx_tc = Wx_tc
        self.Wx_tt = Wx_tt
        self.Wy_t  = Wy_t
        self.Hexp  = Hexp
        self.a0    = a0
        self.Ix    = Ix
        self.Wx    = Wx
        self.rx    = rx
        self.Zx    = Zx
        self.Iy    = Iy
        self.Wy    = Wy
        self.ry    = ry
        self.Cw    = Cw
        self.J     = J
        self.Ag    = Ag
        self.Ix_vz = Ix_vz
        self.Wx_vz = Wx_vz
        self.rx_vz = rx_vz
        self.Zx_vz = Zx_vz
        self.Iy_vz = Iy_vz
        self.Wy_vz = Wy_vz
        self.ry_vz = ry_vz
        self.Cw_vz = Cw_vz
        self.J_vz  = J_vz
        self.Ag_vz = Ag_vz

        self.y_crit      = y_crit
        self.ht_crit     = ht_crit
        self.deffec_crit = deffec_crit
        self.At_crit     = At_crit
        self.yl_crit     = yl_crit
        self.Ix_t_crit   = Ix_t_crit
        self.Iy_t_crit   = Iy_t_crit
        self.rx_t_crit   = rx_t_crit
        self.ry_t_crit   = ry_t_crit
        self.Cw_t_crit   = Cw_t_crit
        self.Zx_t_crit   = Zx_t_crit
        self.J_t_crit    = J_t_crit
        self.Wx_tc_crit  = Wx_tc_crit
        self.Wx_tt_crit  = Wx_tt_crit
        self.Wy_t_crit   = Wy_t_crit
        self.Ix_vz_crit  = Ix_vz_crit
        self.Wx_vz_crit  = Wx_vz_crit
        self.Zx_vz_crit  = Zx_vz_crit   
                     
        self.fya    = dm[nl,0]    
        self.fym    = dm[nl,1]
        self.melast = dm[nl,2]
        self.v      = dm[nl,3]
        self.G      = self.melast/(2*(1+self.v))
                
# Propriedades da seção transversal equivalente (Veríssimo et al. 2012)

        Ie      = 2*(At*y0**2+Ix_t) + (self.tw/24)*(6*Hexp**3+3*Hexp*
                  self.hp**2+8*Hexp**2*self.hp+(2*self.bw/self.p)*
                  (self.hp+Hexp)*(self.hp**2+2*self.hp*Hexp+2*Hexp**2))         # Inércia equivalente   
          
        k1      = (54/(self.tw*y0**2*self.p**2)) * (self.G/self.melast)
        k2      = (0.2*Hexp**3 + 0.375*Hexp*self.hp*(Hexp+0.75*self.hp) + 
                   0.125*self.hp**3)
        k3      = (0.6/(self.tw*y0**2)) * (2.08*Hexp + 1.5*self.hp)
        k4      = (self.p**2 * self.G / (648*self.melast*self.Ix))          
        k5      = (2*self.tw*ya**5/(45*self.Ix**2))         
        
        Ae      = k1 * k2 + k3 + k4 + k5                                        # Área equivalente
        
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.k4 = k4 
        self.k5 = k5
        self.Ie = Ie 
        self.Ae = 1/Ae
        
# Relações da viga

        self.R_LeDG  = int(self.l/self.dg)
        self.R_PeD0  = round((self.p/self.d0),1)
        self.R_DGeD0 = round((self.dg/self.d0),2)
        self.R_D0eD  = round((self.d0/self.d),1)

#----------------------------------------------------------------------------#
# 1.3 - Cálculos Gerais                                                      #
#----------------------------------------------------------------------------#
 
        self.k       = self.dg/self.d                                           # Razão de expansão da viga
        self.np      = self.n-1                                                 # Número de passos da viga
        
#----------------------------------------------------------------------------#
# 1.4 - Seções avaliadas da viga                                             #
#----------------------------------------------------------------------------#
        
# As seções avaliadas para os cálculos dos métodos: 
# As seções são: região dos apoios, centro do montante da alma, 
# centro dos alveolos (aberturas) e plano transversal tangente aos alveolos       
        
        lista        = np.zeros(self.np*4+5)                                    # Vetor das seções avaliadas (4 seções por passo + 5)
        lista[0]     = 0
        lista[1]     = self.bwe                                                 
        
        for i in np.arange(2,4*self.np,4):                                      # Criando as seções que contém em casa passo a partir do primeiro alveolo
            lista[i]   = lista[i-1] + self.a0/2
            lista[i+1] = lista[i]   + self.a0/2
            lista[i+2] = lista[i+1] + self.bw/2
            lista[i+3] = lista[i+2] + self.bw/2
            
            
        lista[self.np*4+2] = self.l-self.bwe-self.a0/2                          # Seção correspondente ao centro do ultimo alveolo
        lista[self.np*4+3] = self.l-self.bwe                                    # Seção correspondente a tangente à direita do ultimo alveolo
        lista[self.np*4+4] = self.l                                             # Seção correspondente ao apoio à direita
        
        self.sections = lista                                         

#----------------------------------------------------------------------------#
# 1.5 - Esforços atuantes                                                    #
#----------------------------------------------------------------------------#
    '''
    Os esforços atuantes pordem ser cargas uniformementes distribuídas ou 
    concentradas no meio do vão da viga

    Também foi considerado apoios do 1º gênero à esquerda e do 2º gênero
    à direita (sistema isostático)
    '''
    def esf(self,q,tipo):                                                       # tipo: caregan=mento uniformemente distribuído (distribuida), carga concentrada (concentrada)
       
        vsd = np.zeros(len(self.sections))
        mf  = np.zeros(len(self.sections))
        
        if tipo  == 'distribuida' :                                             
            Va   = q*self.l/200                                                 # kN/m
        
            for i,x in enumerate(self.sections):                                # enumerate - valores indiciais e valores agregados a cada seção
                vsd[i]  = Va - q * (x/100)                                      # kN
                mf[i]   = Va * (x) - q * (x/100) * (x/2)                        # kN
                            
        if tipo  == 'concentrada' :
            Va       = q/2
            
            for i,x in enumerate(self.sections):
                
                if x < self.l/2:
                    vsd[i] = Va
                    mf[i]  = Va * (x) 
                else:
                    vsd[i]   = Va - q
                    mf[i]    = Va * (x) - q * ((x)-self.l/2)
                    
        self.vsd = abs(vsd)                                                     # Cortante (kN)
        self.mf  = abs(mf)                                                      # Momento fletor (kN.cm)
        self.tipo= tipo                                                         # Tipo de carregamento
        self.q   = q                                                            # Carga atuante no instante de verificação
        
        vsd_ab       = np.zeros(int(self.n/2)+1)                                        
        vsd_ab[0]    = self.vsd[2]                                                   
        vsd_ab[1]    = self.vsd[6]                                                   
        vsd_ab[2:]   = self.vsd[10::4][:int(self.n/2)-1]
        
        vsd_mon       = np.zeros(int(self.n/2)+1)                                        
        vsd_mon[0]    = self.vsd[4]                                                   
        vsd_mon[1]    = self.vsd[8]                                                   
        vsd_mon[2:]   = self.vsd[12::4][:int(self.n/2)-1]
        
        msd_ab       = np.zeros(int(self.n/2)+1)                                        
        msd_ab[0]    = self.mf[2]                                                   
        msd_ab[1]    = self.mf[6]                                                  
        msd_ab[2:]   = self.mf[10::4][:int(self.n/2)-1]
        
        msd_mon       = np.zeros(int(self.n/2)+1)                                        
        msd_mon[0]    = self.mf[4]                                                   
        msd_mon[1]    = self.mf[8]                                                   
        msd_mon[2:]   = self.mf[12::4][:int(self.n/2)-1]
        
        self.vsd_ab   = vsd_ab
        self.vsd_mon  = vsd_mon
        self.msd_ab   = msd_ab
        self.msd_mon  = msd_mon
        
#----------------------------------------------------------------------------#
# 2.0 - Verificação por Veríssimo et al. (2012)                              #
#----------------------------------------------------------------------------# 
            
    def verissimo(self):    
        
# Valores de referência, será substituído ao final das iterações        
        
        ver_FMPS      = 0
        ver_RRS       = 0
        ver_EMAF      = 0
        ver_FMAV      = 0
        ver_FLT       = 0 
        ver_ELS       = 0
        pos_FMPS      = -1                                                               
        pos_RRS       = -1                                                                
        pos_EMAF      = -1                                                                   
        pos_FMAV      = -1                                                             
        pos_FLT       = -1                                                                 
        
#----------------------------------------------------------------------------#
# 2.1 - Formação de Mecanismo Plástico Simples (FMPS)                        #
#----------------------------------------------------------------------------#   
     
        c             = self.y0*self.ya*self.bw*self.At/(2*self.Ix_t)           # Constante que multiplica o esforço cortante (DELESQUES, 1969)
        distesf       = self.vsd*c + self.mf                                    # Distância entre a diferença dos esforços cortante equivalente e momento fletor
       
        distaux       = np.zeros(int(self.n/2)+2)                                        
        distaux[0]    = distesf[2]                                                   
        distaux[1]    = distesf[6]                                                   
        distaux[2:]   = distesf[10::4][:int(self.n/2)]
       
        esfsec        = np.max(distaux)                                         # Maior distância entre os esforços encontrados ao longo da viga
        ind_FMPS      = np.argmax(distaux)                                      # Indice do vetor que indica de maior esforço para falha por FMPS
        pos_FMPS      = self.sections[(ind_FMPS*4)+2]                           # Seção que ocorre falha por FMPS           
        Mpl           = self.fym*self.Zx_vz                                     # Momento de plastificação das seções T's
        Mrd           = Mpl * self.coefseg                                      # Momento resistente (Mpl x Yal) por norma da seção T
        Mesl          = esfsec                                                  # Momento atuante máximo para falha por FMPS
        
        if Mesl > Mrd:                                                          # Condição de falha à FMPS
            ver_FMPS = 'falha'
        
#----------------------------------------------------------------------------#
# 2.2 - Ruptura na Regição da Solda (RRS)                                    #
#----------------------------------------------------------------------------# 
        
        vsdmax        = np.max(self.vsd_mon)                                    # Cortante máximo atuante nos montante de alma na viga
        
        vrk1          = (4/(3*np.sqrt(3))) * (self.bw*self.tw*self.y0*
                         self.fya/self.p)                                       # Esforço cortante resistente característico ao E.M.A.F.
        
        Vrd1          = vrk1 / self.coefseg                                     # Cortante de cálculo ao E.M.A.F.
        
        ind_vsdmax    = np.argmax(vsdmax)                                       # Indice correspondente ao cortante máximo do vetor de cortantes
        pos_vsdmax    = self.sections[(ind_vsdmax+1)*4]                         # Seção correspondente ao cortante máximo
        
        if vsdmax > Vrd1:
            ver_RRS   = 'falha'
            pos_RRS   = pos_vsdmax
            
        self.vrk1 = vrk1            
        
#----------------------------------------------------------------------------#
# 2.3 - Escoamento do Montante da Alma por Flexão (EMAF)                     #
#----------------------------------------------------------------------------#        
        
        eta           = self.p/self.d0                                          # Constante de cálculo
        
        if self.bi != 0 and self.hp >= 0 and self.hp <=                     \
           (self.bw*self.Hexp/(2*self.bi)):                                     # Verificação de viga castelada (se tiver bi = castelada) -> Verificações para cálculo do cortante caracteristico resistente para EMAF
            vrk2      = (8/3) * (self.y0*self.tw*(self.bw*self.Hexp-
                         self.bi*self.hp)*self.fya/(self.hp*self.p))                                                                                                   
        elif self.bi != 0 and self.hp >= (self.bw*self.Hexp/(2*self.bi)): 
            vrk2      = (2/3) * (self.y0*self.tw*self.bw**2*self.fya/
                        (self.hp*self.p))                                                                                           
        else:
            vrk2      = ((self.y0*self.tw*self.fya/(3*eta)) * 
                        ((3*eta-np.sqrt(eta**2+8))/(np.sqrt(4-
                        (eta-np.sqrt(eta**2+8))**2)))) 
        
        Vsd_EMAF       = vrk2 / self.coefseg                                    # Cortante resistente de cálculo referente aos montante da alma entre alveolos
        
        if pos_vsdmax == 0 or pos_vsdmax == self.l:                             # Verificação final do cortente resistente à altura de 0.9R da abertura 
            if self.vsd[4] > Vsd_EMAF:                                          # Verificação na região entre alveolos mais próximo dos apoios
              ver_EMAF   = 'falha'                                                    
              pos_EMAF   = self.sections[4]
        elif vsdmax > Vsd_EMAF:                                                 # Verificação na regiao do cortante máximo fora do apoio
            ver_EMAF   = 'falha'
            pos_EMAF   = self.sections[ind_vsdmax]
            
#----------------------------------------------------------------------------#
# 2.4 - Flambagem do Montante da Alma por Cisalhamento (FMAV)                #
#----------------------------------------------------------------------------# 
                                                                                                                               
        Vcr            = (self.melast*self.tw**3/(1.18*self.y0)) * (1+(1-
                         (2*self.bw/self.p))*((self.y0-0.8*self.Hexp-self.hp)/
                          self.y0))   
        
        if Vcr/vrk2 <= 1:                                                       # Determinação do cortante resistênte à FMAV
            vrk3       = 2/3 * Vcr
        elif Vcr/vrk2 > 1 and Vcr/vrk2 <= 2:
            vrk3       = (vrk2+Vcr)/3
        elif Vcr/vrk2 > 2:
            vrk3       = vrk2
        
        Vsd_FMAV    = vrk3 / self.coefseg                                       # Cortante de cálculo resistente à FMAV na região entre aberturas
 
        for x,i in enumerate(self.vsd_mon[:]):                                  # Looping dos carregamentos para verificação de falha
            if i > Vsd_FMAV:
                ver_FMAV  = 'falha'
                pos_FMAV    = self.sections[(x+1)*4]                            # Posição de falha               
        
#----------------------------------------------------------------------------#
# 2.5 - Flambagem Lateral com Torção (FLT)                                   #
#----------------------------------------------------------------------------#        
     
        mfmax        = np.max(self.mf)                                          # Momento máximo atuante na viga
        tam_sec      = len(self.sections)                                       # Tamanho do vetor das seções analisadas pelas vigas alveolares
        ind_ma       = round(1/4 * tam_sec) - 1                                 # Posição no vetor de seções correspondente à 1/4 da viga
        ma           = self.mf[ind_ma]                                          # Momento atuante na posição 1/4 do comprimento da viga
        ind_mb       = round(2/4 * tam_sec)                                     # Posição no vetor das seções correspondente ao meio da viga 
        mb           = self.mf[ind_mb]                                          # Momento atuante no meio da viga
        ind_mc       = round(3/4 * tam_sec)                                     # Posição no vetor das seções corresponte à 3/4 da viga 
        mc           = self.mf[ind_mc]                                          # Momento atuante à 3/4 da viga
        cb           = 12.5*mfmax / (2.5*mfmax + 3*ma + 4*mb + 3*mc)            # Constante de cálculo (Usado para determinar o momento característico)
        beta1        = 0.7*self.fym*self.Wx_vz/(self.melast*self.J_vz)          # Constante de cálculo
        Lp           = (1.76*self.ry_vz*np.sqrt(self.melast/self.fym))          # Comprimento de plastificação                                                                
        Lrcor        = (1.66*np.sqrt(self.Iy_vz*self.J_vz)/
                       (self.J_vz*beta1))* (np.sqrt(1+np.sqrt(1+
                       (27*self.Cw*beta1**2/self.Iy_vz))))                      # Comprimento de escoamento
        Mrcor        = (0.31*self.melast/Lrcor**2) * (np.sqrt(self.Iy_vz*
                       (1000*self.Cw + 39*self.J_vz*self.lb**2)))               # Momento de início de escoamento
        
        if self.lb > Lrcor :                                                    # Determinação do momento resistente caracteristico   
            Mrk      = ((cb*np.pi**2*self.melast*self.Iy_vz/self.lb**2) *
                       (np.sqrt((self.Cw/self.Iy_vz) * 
                       (1+0.39*self.J_vz*self.lb**2/self.Cw))))
        elif Lp < self.lb and self.lb <= Lrcor:
            Mrk      = (cb * (0.9 * Mpl - (0.9 * Mpl - Mrcor) * 
                       ((self.lb - Lp)/(Lrcor - Lp))))
        else:
            Mrk      = 0.9 * Mpl
        
        ind_FLT      = np.argmax(self.mf)                                       # Indice correspondente ao momento máximo do vetor de cortantes  
        pos_FLT      = self.sections[ind_FLT]                                   # Posição que ocorre o momento máximo
        
        Msd_FLT      = Mrk / self.coefseg                                       # Momento resistente à FLT
        
        if mfmax > Msd_FLT:                                                     # Verificação de falha
            ver_FLT  = 'falha'
            
#----------------------------------------------------------------------------#
# 2.6 - Verificação de Flecha Máxima (ELS)                                   #
#----------------------------------------------------------------------------#
        
# Flecha máxima foi baseada na formulação proposta por TIMOSHENKO (1966).
# Os deslocamentos máximos permitidos foram como base a NBR 8800:2008

        if self.tipo == 'distribuida':
            w         = self.q/100
            Um        = (5/384) * w* self.l**4 / self.melast / self.Ie
            Uad       = w * self.l**2 / 8 / self.G / self.Ae
            Uv        = 0
        elif self.tipo == 'concentrada': 
            alfa      = (self.Ae/(self.tw*self.Ie)) * ((self.bf*self.dg**2/8) 
                        - ((self.dg-2*self.tf)**2/8) * (self.bf-self.tw))
            Uad       = (alfa / self.Ae / self.G) * (self.q * self.l / 4)
            Uv        = (self.q * self.l**3 / (48 * self.melast * self.Ie))
            Um        = 0
                    
        Umax          = (Um + Uv + Uad)
                
        if self.sist  == 'piso':
            Ulim      = self.l / 350
        else:
            Ulim      = self.l / 250
                
        if Umax > Ulim:
            ver_ELS  = 'flecha excedida'
        else:
            ver_ELS  = 'OK'      
 
#----------------------------------------------------------------------------#
# 2.7 - Atribuições do método Veríssimo ao objeti viga alveolar              #
#----------------------------------------------------------------------------#
                       
        self.pos_FMPS = pos_FMPS
        self.esfsec   = esfsec
        self.pos_RRS  = pos_RRS
        self.pos_EMAF  = pos_EMAF
        self.pos_FMAV = pos_FMAV
        self.pos_FLT  = pos_FLT
        self.ver_ELS  = ver_ELS
        self.flecha   = Umax
        self.flechalim= Ulim
        
        return ver_RRS,ver_FMPS,ver_EMAF,ver_FMAV,ver_FLT,ver_ELS
    
#----------------------------------------------------------------------------#
# 3.0 - Verificação por AISC 31 (2016)                                       #
#----------------------------------------------------------------------------#     

    def AISC(self):
        
#----------------------------------------------------------------------------#
# 3.1 - Formação de Mecanismo de Vierendeel (FMV)                            #
#----------------------------------------------------------------------------#        
        
# Esforços Requeridos
        
        if self.bi != 0:
            Pr       = self.msd_ab / self.deffec                                    # Solicitação axial na seção
        else:
            Pr       = self.msd_ab / self.deffec_crit                               # Solicitação axial na seção (celular)
            
        if self.bi != 0:
            Mvr  = self.vsd_ab * (self.At / self.Anet) * (self.bw / 2)             # Momento de Vierendeel solicitado na seção
        else:
            Mvr  = self.vsd_ab * (1 / 2) * (self.d0 / 4)
        
# Dados para verificação de falha por Vierendeel
        
        kx       = 0.65                                                         # Coeficiente (X)
        ky       = 1.00                                                         # Coeficiente (Y) 
        kz       = 1.00                                                         # Coeficiente (Z)

        if self.bi != 0:                                                        # Verificação se a seção é castelar ou celular
            lld  = self.bw                                                      # Comprimento lateralmente destravado na região dos T's (Castelada)
        else:
            lld  = self.d0/2                                                    # Comprimento lateralmente destravado na região dos T's (Celular)
            
        lcx      = kx * lld                                                     # Coeficiente de comprimento destravado em relação a X
        lcy      = ky * lld                                                     # Coeficiente de comprimento destravado em relação a Y
        lcz      = kz * lld                                                     # Coeficiente de comprimento destravado em relação a Z
        
        rmin_t   = min(self.ry_t,self.rx_t)                                     # Menor raio de giração da seção
        
# Verificação da resistência a compressão para flambagem por flexão (1)
        
        ax        = lcx/self.rx_t
        by        = lcy/self.ry_t
        
        const1   = max(ax,by) / rmin_t                                          # Constante da razão de limite para cálculo da força crítica
        
        Fe1       = np.pi**2 * self.melast / const1**2                          # Força elástica para resistencia a compressão
        
        if const1 <= (4.71*np.sqrt(self.melast/self.fya)) or                \
           (self.fya/Fe1) <= 2.25:                                              # Deterimação da força crítica para este modo de falha
            Fcr1  = (0.658**(self.fya/Fe1)) * self.fya
        elif const1 > (4.71*np.sqrt(self.melast/self.fya)) or               \
             (self.fya/Fe1) > 2.25:
            Fcr1  = 0.877 * Fe1
        
        Pn1      = (Fcr1 * self.At_crit) / self.coefseg                         # Resistência à compressão para flambagem à flexão dos T's 
        
# Verificação da resistência a compressão para flambagem por flexo-torção (2)
        
        if self.bi != 0:
            xo         = 0      
            yo         = (self.ht - self.yl) - (self.tf/2)     
        
            ro2        = (xo**2 + yo**2 + ((self.Ix_t+
                                 self.Iy_t)/self.At))                           # Raio polar de giração da seção T 
            Fey        = np.pi**2 * self.melast / (lcy/self.ry_t)**2            # Força elástica na direção Y
            
            Fey_ksi    = Fey*1.4503773772954
            
            Fez_ksi    = (((np.pi**2*self.melast*1.4503773772954/
                         ((lcz/2.54)**2)) + self.G*1.4503773772954 * 
                           (self.J_vz*0.0240251480/2)) / 
                          (self.At*ro2*0.155**2))                               # Força elástica na direção X
            
            Fez        = Fez_ksi *  0.6894757
                
            H          = 1 - (xo**2 + yo**2) / ro2                              # Constante de cálculo
            Fe2_ksi    = ((Fey_ksi+Fez_ksi)/(2*H)) * ( 1 - np.sqrt(1 - 
                         (4*Fey_ksi*Fez_ksi*H)/(Fey_ksi+Fez_ksi)**2))
        
            Fe2        = Fe2_ksi *  0.6894757
        else:
            xo         = 0                  
            yo         = (self.ht_crit - self.yl_crit) - (self.tf/2)

            ro2        = (xo**2 + yo**2 + ((self.Ix_t_crit+
                                 self.Iy_t)/self.At_crit))                      # Raio polar de giração da seção T 
            Fey        = np.pi**2 * self.melast / (lcy/self.ry_t_crit)**2       # Força elástica na direção Y
            
            Fey_ksi    = Fey*1.4503773772954
            
            Fez_ksi    = (((np.pi**2*self.melast*1.4503773772954/
                         ((lcz/2.54)**2)) + self.G*1.4503773772954 * 
                           (self.J_vz*0.0240251480/2)) / 
                          (self.At_crit*ro2*0.155**2))                          # Força elástica na direção X
            
            Fez        = Fez_ksi *  0.6894757
                
            H          = 1 - (xo**2 + yo**2) / ro2                              # Constante de cálculo
            Fe2_ksi    = ((Fey_ksi+Fez_ksi)/(2*H)) * ( 1 - np.sqrt(1 - 
                         (4*Fey_ksi*Fez_ksi*H)/(Fey_ksi+Fez_ksi)**2))
        
            Fe2        = Fe2_ksi *  0.6894757
        
        if const1 <= (4.71*np.sqrt(self.melast/self.fya)) or                \
           (self.fya/Fe2) <= 2.25:                                              # Determinação da força crítica para compressão à flambagem por flexo-torção
            Fcr2  = (0.658**(self.fya/Fe2)) * self.fya
        elif const1 > (4.71*np.sqrt(self.melast/self.fya)) or               \
           (self.fya/Fe2) > 2.25:
            Fcr2  = 0.877 * Fe2
        
        Pn2      = (Fcr2 * self.At_crit) / self.coefseg                         # Resistência a compressão para flambagem por flexo- compressão        
        
# Verificação à resistência a tração  
        
        Pn3      = self.fya * self.At_crit / self.coefseg                       # Resistencia à Tração
        
# Resistência normal axial dos T's
        
        Pn       = min(Pn1,Pn2,Pn3)                                             # Resistencia mínima axial para falha por FMV  
        
# Resistência nominal a flexo plastificação (1) 
        
        My       = self.fya * self.Wx_tc_crit                                   # Momento de início de plastificação
        Mp       = My                                    
        Mn1      = Mp 
        
# Resistência a flambagem Lateral-Torsão (2)
        
        Lb_B     = lld                                                          # Comprimento destravado na região do T 
        Lp_B     = 1.76*self.ry_t*np.sqrt(self.melast/self.fya)                 # Comprimento máximo referente à falha por plastificação do T
        Lr_B     = 1.95*(self.melast/self.fya) * (np.sqrt(self.Iy_t*self.J_t)/
                   self.Wx_tc_crit) * np.sqrt(2.36*(self.fya/self.melast)*
                   (self.ht_crit*self.Wx_tc_crit)/self.J_t_crit + 1)            # Comprimento máximo referente à falha por Flambagem inelástica do T
        B        = (2.3*self.ht_crit/Lb_B)*np.sqrt(self.Iy_t_crit/
                                                   self.J_t_crit)               # Constante de Cálculo - Em compressão
        Mcr2     = ((1.95*self.melast/Lb_B) * np.sqrt(self.Iy_t_crit*
                    self.J_t_crit) * (B + np.sqrt(1 + B**2)))                   # Momento crítico à FLT
        
        if Mcr2 > My:                                                           # Limite de falor para o Momento crítico
            Mcr2 = My
        
        if Lb_B <= Lp_B:                                                        # Verificação do momento resistente a flambagem lateral-torção
            Mn2  = 9e99                                                         # Valor alto para inticar que: Não se aplica 
        elif Lp_B < Lb_B and Lb_B <= Lr_B:
            Mn2  = Mp - (Mp-My) - ((Lb_B - Lp_B)/(Lr_B - Lp_B))
        elif Lb_B > Lr_B:
            Mn2  = Mcr2
        
# Resistência de flambagem local nos T's (3)
       
        lamb      = self.bf / (2*self.tf)                                       # Índice de esbeltez da Mesa
        lamb_p    = 0.38 * np.sqrt(self.melast/self.fya)                        # Índice de esbeltez referente à falha por plastificação da mesa
        lamb_r    = 1.00 * np.sqrt(self.melast/self.fya)                        # Índice de esbetez referente à falha por flambagem Inelástica da Mesa
        
        if lamb <= lamb_p:                                                      # Verificação do momento resistente a flambagem local dos T's
            Mn3   = 9e99                                                        # Valor alto para inticar que: Não se aplica 
            comp  = 'compacta'
        elif lamb > lamb_p and lamb <= lamb_r:
            Mn3   = (Mp - (Mp - 0.75*self.fya*self.Wx_tc_crit) * 
                    ((lamb - lamb_p)/(lamb_r - lamb_p)))
            comp  = 'não compacta'
            if Mn3 > 1.6*My:
                Mn3  = 1.6*My
        elif lamb > lamb_r:
            Mn3   = 0.75*self.melast*self.Wx_tc_crit/lamb**2
            comp  = 'slender'
            
# Flambagem Local de Hastes em T na Flexo-compressão (4)
        
        dt        = (self.dg - self.d0) / 2                                     # Distância do topo da viga até a abertura
        
        if dt/self.tw <= 0.84*np.sqrt(self.melast/self.fya):                    # Determinação da força crítica 
            Fcr4  = self.fya
        elif 0.84*np.sqrt(self.melast/self.fya) < dt/self.tw and            \
             dt/self.tw <= 1.52*np.sqrt(self.melast/self.fya):
            Fcr4  = (1.43 - 0.515 * (dt/self.tw) * 
                     np.sqrt(self.fya/self.melast)) * self.fya
        elif dt/self.tw > 1.52*np.sqrt(self.melast/self.fya):    
            Fcr4  = 1.52 * self.melast / (dt/self.tw)**2
            
        Mn4       = Fcr4 * self.Wx_tc_crit                                      # Momento de resistência à flambagem local no T na flexo-compressão
        
# Resistência normal flexão dos T's 
        
        Mn        = min(Mn1,Mn2,Mn3,Mn4)/ self.coefseg                          # Momento resistivo final para FMV nos T's               
        
# Verificação de Resistência Axial 
        
        AISC_FMVa       = 0                                                     # Valor de referência, será substituído posteriormente
        
        for x,i in enumerate(Pr[:]):                                            # Looping entre as seções das vigas 
            if i > Pn:                                                          # Verificação da solicitação X resistencia  
                AISC_FMVa       = 'falha'
                B_pos_FMVa      = self.sections[(x*2)+1]                        # Definição da posição de falha 
                self.B_pos_FMVa = B_pos_FMVa
         
# Verificação de Resistência Axial 
        
        AISC_FMVm       = 0                                                     # Valor de referência, será substituído posteriormente
        
        for x,i in enumerate(Mvr[:]):                                           # Looping entre as seções das vigas 
            if i > Mn:                                                          # Verificação da solicitação X resistencia 
                AISC_FMVm       = 'falha'
                B_pos_FMVm      = self.sections[(x*2)+1]                        # Definição da posição de falha 
                self.B_pos_FMVm = B_pos_FMVm
                
 # Verificação da iteração entre a flexão e força axial nos T's   
        
        AISC_FMV  = 0                                                           # Valor de referência, será substituído posteriormente
        B_pos_FMV = -1                                                          # Valor de referência, será substituído posteriormente 
        const2    = np.zeros(len(Pr))                                           # Valor de referência, será substituído posteriormente
        
        for i in range(len(Pr)):                                                # Looping entre as seções das vigas
            if Pr[i] / Pn >= 0.2:                                               # Verificação da iteração entre os esforços
                const2[i]  = (Pr[i] / Pn) + (8/9) * (Mvr[i] / Mn)               # Relação da iteração dos esforços
                if const2[i] > 1.0:                                             # Verificação de falha
                    AISC_FMV       = 'falha'
                    B_pos_FMV      = self.sections[(i*2)+1]                     # Posição de falha
                    self.B_pos_FMV = B_pos_FMV
            elif Pr[i] / Pn < 0.2:    
                const2[i]  = (Pr[i] / (2*Pn)) + (Mvr[i] / Mn)                   # Relação da iteração dos esforços
                if const2[i] > 1.0:                                             # Verificação de falha
                    AISC_FMV       = 'falha'
                    B_pos_FMV      = self.sections[(i*2)+1]                     # Posição de falha
                    self.B_pos_FMV = B_pos_FMV 
        
#----------------------------------------------------------------------------#
# 3.2 - Ruptura na Regição da Solda (RRS)                                    #
#----------------------------------------------------------------------------#
 
        Vrh_ap   = self.vsd[0] * (self.p / (self.dg - 2*self.yl))

        Vn      = 0.6*self.fya*self.bw*self.tw                                  # Resistência nominal de corte horizontal na alma 
        Vn_ap   = 0.6*self.fya*self.bwe*self.tw                                 # Resistência nominal de corte horizontal na alma (região do apoio) 

        mdiff = np.sqrt((self.msd_ab[1:] - self.msd_ab[:-1])**2)                # Diferença dos momentos atuantes à direita e à esquerda do montante da alma
        
        Vrh   = mdiff / self.deffec                                             # Cortante horizontala atuante na solda 
                
        AISC_RRS = 'não falha'                                                  # Atribuição inicial
        if Vrh_ap > Vn_ap:                                                      # Verificação de falha na região próxima ao apoio
            AISC_RRS   = 'falha'
            pos_RRS   = 0                                                       # Posição de falha no apoio
            self.pos_RRS = pos_RRS    
        else:                                                                   # Verificação de falha nas regigçoes entre alveolos
            for x,i in enumerate(Vrh[:]):                                       # Looping dos carregamentos para verificação de falha
                if i > Vn:
                    AISC_RRS   = 'falha'
                    B_pos_RRS    = self.sections[(x+1)*4]                       # Posição de falha
                    self.B_pos_RRS    = B_pos_RRS                        
         
#----------------------------------------------------------------------------#
# 3.3 - Flambagem do Montante da Alma por Cisalhamento (FMA)                 #
#----------------------------------------------------------------------------#
                    
        Mpbm     = 0.25 * self.tw * (self.bw + 2*self.bi)**2 * self.fya         # Momento de flexão plástico 
        Mpbm_ap  = 0.25 * self.tw * (self.bwe + 2*self.bi)**2 * self.fya        # Momento de flexão plástico (região de apoio)
        Mebm     = self.tw*(self.p - self.d0 + 0.564*self.d0)**2 * self.fya/6   # Momento de flexão elástico 
        Mebm_ap  = self.tw*((self.p-self.bw+self.bwe) - self.d0 +
                   0.564*self.d0)**2 * self.fya / 6                             # Momento de flexão elástico (região de apoio)
        r1       = self.bw/self.tw                                              # Relação de largura do alveolo e espessura da alma   
        B_cond     = 'Aplica'                                                   # Condição de aplicabilidade do método B para vigas casteladas    
        
# Solicitações nas viga 

        if self.bi != 0:                                                        # Verificação de viga castelada 
            Mrh    = Vrh * self.ht 
            Mrh_ap = Vrh_ap * self.ht                                           # Momento solicitante para vigas castelada
        elif self.bi == 0:                                                      # Verificação de viga celular
            Mrh = 0.9 * self.d0 * Vrh / 2                                       # Momento solicitante para viga celular
            Mrh_ap = 0.9 * self.d0 * Vrh_ap / 2

# Verificação em vigas casteladas
        
        # Altura do lado inclinado da abertura castelada
        
        if self.hp != 0:                                                        # Verificação se há chapa espansiva 
            h_inc= self.h0/2 - self.hp/2                                        # Altura da parte inclinada se houver chapa espansiva
        else:
            h_inc= self.h0/2                                                    # Altura da parte inclinada se não houver chapa espansiva
        
        if self.bi != 0:                                                        # Verificação de viga castelada
            teta     = np.degrees(np.arctan(h_inc/self.bi))                     # Angulo formado entre a parte inclinada e a linha horizontal
        else:
            teta     = 'Não possui'                                             # Não há ângulo para vigas celulares                             
        
        # Ajuste de relação entre lagura das aberturas e espessura da alma
        
        if r1 >= 7 and r1 <= 13:                                                # r_aj é um ajuste feito para a relação r1, assim como o guia americano aproxima a angulação, foi feita essa aproximação para que não seja tão específico e deterministico
            r_aj = 10
        elif r1 >= 17 and r1 <= 23:
            r_aj = 20
        elif r1 >= 27 and r1 <= 33:
            r_aj = 30
        else:
            B_cond = 'Não aplica'                                               # Indicando que o método B não atende as condições de aplicabilidade para vigas casteladas 
            r_aj   = 0                                                          # Assumindo valor fora dos limites
        
        if self.bi != 0:                                                        # Verificação de viga castelada
            
            if teta >= 43 and teta <=47:                                        # Verificação do angulo entres os valores limites para formulação de teta = 45º
                
                if r_aj == 10:                                                  
                    Mocr    = Mpbm * (0.351 - 0.051 * (h_inc*2/self.bw) +       # Relação de largura do montante da alma e espessura = 10 
                              0.0026 * (h_inc*2/self.bw)**2)                    # Momento crítico para verificação de FMA em vigas casteladas
                    Mocr_ap = Mpbm_ap * (0.351 - 0.051 * (h_inc*2/self.bwe) +
                              0.0026 * (h_inc*2/self.bwe)**2)                   # Momento crítico para verificação de FMA em vigas casteladas (Região dos apoios)
                elif r_aj == 20:                                                # Relação de largura do montante da alma e espessura = 20 
                    Mocr    = Mpbm * (3.276 - 1.208 * (h_inc*2/self.bw) + 
                              0.154 * (h_inc*2/self.bw)**2 -
                              0.0067 * (h_inc*2/self.bw)**3)                    # Momento crítico para verificação de FMA em vigas casteladas
                    Mocr_ap = Mpbm_ap * (3.276 - 1.208 * (h_inc*2/self.bwe) + 
                              0.154 * (h_inc*2/self.bwe)**2 - 
                              0.0067 * (h_inc*2/self.bwe)**3)                   # Momento crítico para verificação de FMA em vigas casteladas (Região dos apoios)
                elif r_aj == 30:                                                # Relação de largura do montante da alma e espessura = 30 
                    Mocr    = Mpbm * (0.952 - 0.30 * (h_inc*2/self.bw) + 
                                      0.0319 * (h_inc*2/self.bw)**2 - 
                                      0.0011 * (h_inc*2/self.bw)**3)            # Momento crítico para verificação de FMA em vigas casteladas
                    Mocr_ap = Mpbm_ap * (0.952 - 0.30 * (h_inc*2/self.bwe) + 
                                         0.0319 * (h_inc*2/self.bwe)**2 - 
                                         0.0011 * (h_inc*2/self.bwe)**3)        # Momento crítico para verificação de FMA em vigas casteladas (Região dos apoios)
            
            elif  teta >= 58 and teta <=62:                                     # Verificação do angulo entres os valores limites para formulação de teta = 60º
                
                if r_aj == 10:                                                  # Relação de largura do montante da alma e espessura = 10 
                    Mocr    = Mpbm * (0.587 * (0.917)**(2*h_inc/self.bw))       # Momento crítico para verificação de FMA em vigas casteladas
                    Mocr_ap = Mpbm_ap * (0.587 * (0.917)**(2*h_inc/self.bwe))   # Momento crítico para verificação de FMA em vigas casteladas (Região de apoios)
                elif r_aj == 20:                                                # Relação de largura do montante da alma e espessura = 20 
                    Mocr    = Mpbm * (1.960 * (0.699)**(2*h_inc/self.bw))       # Momento crítico para verificação de FMA em vigas casteladas  
                    Mocr_ap = Mpbm_ap * (1.960 * (0.699)**(2*h_inc/self.bwe))   # Momento crítico para verificação de FMA em vigas casteladas (Região de apoios)
                elif r_aj == 30:                                                # Relação de largura do montante da alma e espessura = 30 
                    Mocr    = Mpbm * (2.55 * (0.574)**(2*h_inc/self.bw))        # Momento crítico para verificação de FMA em vigas casteladas
                    Mocr_ap = Mpbm_ap * (2.55 * (0.574)**(2*h_inc/self.bwe))    # Momento crítico para verificação de FMA em vigas casteladas (Região de apoios)
            
            else: 
                
                B_cond = 'Não aplica'                                           # Verificação não aplicabilidade do método B para vigas casteladas
            
# Verificação em vigas celulares

        C1       = (5.097 + 0.1464 * (self.d0/self.tw) - 
                    0.00174  * (self.d0/self.tw)**2)                               
        C2       = (1.441 + 0.0625 * (self.d0/self.tw) - 
                    0.000683 * (self.d0/self.tw)**2)                            # Constante de calculo 2 para vigas celulares 
        C3       = (3.645 + 0.0853 * (self.d0/self.tw) -
                    0.00108  * (self.d0/self.tw)**2)                            # Constante de calculo 3 para vigas celulares 
         
        Mallow   = Mebm * (C1 * (self.p/self.d0) - 
                   C2 * (self.p/self.d0)**2 - C3)                               # Momento resistivo para vigas celulares
        Mallow_ap= Mebm_ap * (C1 * (self.p/self.d0) - 
                   C2 * (self.p/self.d0)**2 - C3)                               # Momento resistivo para vigas celulares (Região dos apoios) 
     
        if self.bi == 0:                                                        # Verificação de vigas celulares
            B_cond  = 'Aplica'                                                  # Verificação de aplicabilidade assumindo ok para vigas celulares
       
# Verificação final

        if self.bi != 0:                                                        # Verificação de viga castelada
            Mn2    = Mocr / self.coefseg                                        # Momento de resistencia de projeto para falha por FMA
            Mn2_ap = Mocr_ap / self.coefseg                                     # Momento de resistencia de projeto para falha por FMA (Região de apoio)
        elif self.bi == 0:                                                      # Verificação de viga celular
            Mn2    = Mallow / self.coefseg                                      # Momento de resistencia de projeto para falha por FMA
            Mn2_ap = Mallow_ap/self.coefseg                                     # Momento de resistencia de projeto para falha por FMA (Região de apoio)  
            
        AISC_FMA = 'não falha'                                                  # Valor de referência para caso não entre no lopping
    
        for x,i in enumerate(Mrh[:]):
            if i > Mn2:
                AISC_FMA      = 'falha'
                B_pos_FMA     = self.sections[(x+1)*4]                          # Determinação da posição de falha associando o vetor de solicitações com o vetor das seções
                self.B_pos_FMA= B_pos_FMA    
#----------------------------------------------------------------------------#
# 3.4 - Flambagem do Montante da Alma por Compressão (FMAC)                  #
#----------------------------------------------------------------------------#

        kv_gross  = 5.34                                                        # Constante de cálculo para seção de viga de alma cheia
        kv_net    = 1.20                                                        # Constante de cálculo para seção de viga com abertura
        h_gross   = self.dg - 2 * self.tf                                       # Altura da alma 
        dt        = self.ht - self.tf                                           # Metade da altura da alma na região de abertura
        r_gross   = h_gross / self.tw                                           # Relação de altura da alma com a espessura alma na região do montante de alma cheia 
        r_net     = dt / self.tw                                                # Relação de altura da alma com a espessura alma na região do alveolo  
        
# Determinação dos coeficientes de cálculo        
        
        if r_gross <= 1.1 * np.sqrt(kv_gross*self.melast/self.fya):
            Cv1   = 1.0
        elif r_gross > 1.1 * np.sqrt(kv_gross*self.melast/self.fya):
            Cv1   = (1.1 * np.sqrt(kv_gross*self.melast/self.fya)) / r_gross
         
        if r_net <= 1.1 * np.sqrt(kv_net*self.melast/self.fya): 
            Cv2   = 1.0
        elif r_net > 1.1 * np.sqrt(kv_net*self.melast/self.fya) \
             and r_net <= 1.37 * np.sqrt(kv_net*self.melast/self.fya):
            Cv2   = (1.1 * np.sqrt(kv_net*self.melast/self.fya)) / r_net
        elif r_net > 1.37 * np.sqrt(kv_net*self.melast/self.fya):
            Cv2   = 1.5 * kv_net * self.melast / (r_net**2 * self.fya)
                         
        Vn_gross  = 0.6 * self.fya * self.dg * self.tw * Cv1                    # Resistência nominal de corte na regiao do montante de alma cheia
        Vn_net    = 0.6 * self.fya * (2*dt) * self.tw * Cv2                     # Resistência nominal de corte na região das aberturas
        
# Solicitações 

        vsd_gross = np.zeros(int(self.n)+1)                                     # Criação do vetor das solicitações na região do montante de alma cheia
        vsd_gross[0]   = self.vsd[1]                                            # Primeira posição sendo a seção transversal da viga tangente a primeira abertura
        vsd_gross[1:-1]= self.vsd[4::4][:(int(self.n)-1)]                       # Demais seções do montante de alma cheia
        vsd_gross[-1]  = self.vsd[-2]                                           # Ultima posição sendo a seção transversal da tangente a ultima abertura

        vsd_net = np.zeros(int(self.n))                                         # Criação dp vetor das solicitações na região onde há abertura
        vsd_net[0] = self.vsd[2]                                                # Primeira seção sendo correspondente à primeira abertura
        vsd_net[1:]= self.vsd[6::4][:(int(self.n)-1)]                           # Demais seções transversais das demais aberturas

# Verificação final

        AISC_FMAC     = 'não falha'                                             # Valor de referência, depois será substituído
        AISC_FMACnet  = 'não falha'                                             # Valor de referência, depois será substituíd
        B_pos_FMAC    = -1                                                      # Valor de referência, depois será substituíd
        B_pos_FMACnet = -1                                                      # Valor de referência, depois será substituíd
        
        for x,i in enumerate(vsd_gross):                                        # Lopping para verificação final 
            if i > Vn_gross:                                                    # Comparação do esforço cortante na seção em avaliação com a resistênte 
                if x == 0:                                                      # Posção referente a como sendo próxima ao apoio da esquerda
                    AISC_FMAC       = 'falha'                                   
                    B_pos_FMAC      = self.sections[1] 
                    break
                elif x == len(vsd_gross)-1:                                     # Posição referente como sendo próxima ao apoio da direita
                    AISC_FMAC       = 'falha'
                    B_pos_FMAC      = self.sections[-2]
                    break
                else:                                                           # Demais seções ao longo da viga correspondente aos montante de alma cheia
                    AISC_FMAC       = 'falha'
                    B_pos_FMAC      = self.sections[(4*x)]
                    break
        self.B_pos_FMAC             = B_pos_FMAC
        
        for x,i in enumerate(vsd_net):                                          # Lopping para verificação final da região onde há abertura 
            if i > Vn_net:
                if x == 0:
                    AISC_FMACnet    = 'falha'
                    B_pos_FMACnet   = self.sections[2]
                    break
                else:    
                    AISC_FMACnet    = 'falha'
                    B_pos_FMACnet   = self.sections[(4*x)]
                    break
        self.B_pos_FMACnet          = B_pos_FMACnet

#----------------------------------------------------------------------------#
# 3.5 - Flambagem Lateral por Torção (FLT)                                   #
#----------------------------------------------------------------------------#

        mfmax        = np.max(self.mf)                                          # Momento máximo atuante na viga                                                                                      
        tam_sec      = len(self.sections)                                       # Número de seções avalidadas ao longo da viga                                                              
        ind_ma       = round(1/4 * tam_sec) - 1                                 # Indice do vertor das seções correspondente à posição 1/4 do comprimento                                                         
        ma           = self.mf[ind_ma]                                          # Momento atuante na posição 1/4 do comprimento da viga
        ind_mb       = round(2/4 * tam_sec)                                     # Posição no vetor das seções correspondente ao meio da viga 
        mb           = self.mf[ind_mb]                                          # Momento atuante no meio da viga
        ind_mc       = round(3/4 * tam_sec)                                     # Posição no vetor das seções corresponte à 3/4 da viga 
        mc           = self.mf[ind_mc]                                          # Momento atuante à 3/4 da viga
        cb           = 12.5*mfmax / (2.5*mfmax + 3*ma + 4*mb + 3*mc)            # Fator de modificação para falha lateral com torção   
        AISC_FLT     = 'não falha'                                              # Valor de referência, depois será alterado
        B_pos_FLT    = -1                                                       # Valor de referência, depois será alterado
        
        if comp == 'compacta':                                                  # Verificações em situação de viga compacta
            Mp_FLT   = self.fym * self.Zx                                       # Momento plástico
            c        = 1                                                        # Constante de dupla simetria
            rts      = np.sqrt(np.sqrt(self.Iy*self.Cw) / self.Wx)              # Raio de giração
            Lp       = 1.76 * self.ry * np.sqrt(self.melast/self.fya)           # Comprimento máximo referente à falha por plastificação 
            Lr       = 1.95 * rts * (self.melast/(0.7*self.fya)) *          \
                       np.sqrt((self.J*c/self.Wx/(self.dg-2*self.tf)) +
                       np.sqrt((self.J*c/self.Wx/(self.dg-
                       2*self.tf))**2 + 6.76 * (0.7*self.fya/self.melast)**2))  # Comprimento máximo referente à falha por flambagem
            Fcr      = ((cb * np.pi**2 * self.melast /(self.lb/rts)**2) *     
                       np.sqrt(1 + 0.078 * (self.J*c/self.Wx/               
                       (self.dg-2*self.tf)) * ((self.lb/rts)**2)))              # Força crítica
            
            if self.lb <= Lp:                                                   # Determinação do momneto nominal resistivo para FLT quando a seção é compacta                                   
                Mn_FLT = 9e99                                                   # Valor alto para indicar que não se aplica
            elif Lp < self.lb <= Lr:
                Mn_FLT = cb * (Mp_FLT - (Mp_FLT - 0.7*self.fym * self.Wx) * 
                         ((self.lb - Lp)/(Lr - Lp)))
            elif self.lb > Lr:
                Mn_FLT = Fcr * self.Wx
                
            for x,i in enumerate(self.mf):                                      # Verificação final de falha
                if i > Mn_FLT:
                    AISC_FLT    = 'falha'
                    B_pos_FLT   = self.sections[x]                              # Posição de falha
                    # break

        elif comp == 'não compacta':                                            # Verificações em situação de viga não compacta
             Mp_FLT   = self.fym * self.Zx                                      # Momento de plastificação
             Mn_FLT   = (Mp_FLT - (Mp_FLT - 0.7*self.fym * self.Wx) * 
                        ((lamb - lamb_p)/(lamb_r - lamb_p)))                    # Determinação do momneto nominal resistivo para FLT quando a seção é não compacta
             
             for x,i in enumerate(self.mf):                                     # Verificação final de falha
                 if i > Mn_FLT:
                     AISC_FLT    = 'falha'
                     B_pos_FLT   = self.sections[x]                             # Posição de falha
        
        elif comp == 'slender':                                                 # Verificações em situação de viga é esbelta
            kc       = 4 / np.sqrt((self.dg-2*self.tf))                         # Constante de cálculo
            Mp_FLT   = self.fym * self.Zx                                       # Momento plástico
            Mn_FLT   = 0.9 * self.melast * kc * self.Wx / lamb**2               # Determinação do momneto nominal resistivo para FLT quando a seção é esbelta
            
            for x,i in enumerate(self.mf):                                      # Verificação de falha à FLT quando a seção é esbelta
                if i > Mn_FLT:
                    AISC_FLT    = 'falha'      
        
#----------------------------------------------------------------------------#
# 3.6 - Verificação de Flecha Máxima (ELS)                                   #
#----------------------------------------------------------------------------#

        I_AISC        = 0.9 * self.Ix_vz                                        # Simplificação feita pelo guia americano AISC 31
        A_AISC        = 0.9 * self.Ag                                           # Simplificação feita pelo guia americano AISC 31
        
        if self.tipo == 'distribuida':
            w         = self.q/100                                              # Transformando de kN/m para kN/cm
            Um        = (5/384) * w* self.l**4 / self.melast / I_AISC           # Deslocamento devido ao carregamento uniformemente distribuído
            Uv        = 0                                                       
            Uad       = w * self.l**2 / 8 / self.G / self.Ae                    # Deslocamento adicional devido ao cortante pelo carregamento uniformemente distribuído 
        elif self.tipo == 'concentrada': 
            alfa      = (A_AISC/(self.tw*I_AISC)) * ((self.bf*self.dg**2/8) -   # Constante de cálculo (TIMOSHENKO, 1966)
                        ((self.dg-2*self.tf)**2/8) * (self.bf-self.tw))
            Uad       = (alfa / A_AISC / self.G) * (self.q * self.l / 4)        # Deslocamento adicional devido ao cortante por cargas concentradas (TIMOSHENKO, 1966)
            Uv        = (self.q * self.l**3 / (48 * self.melast * I_AISC))      # Deslocamento devido ao momento provocado por cargas concentradas (TIMOSHENKO, 1966)
            Um        = 0
            
        Umax          = (Um) #+ Uv + Uad)                                         # Deslocmaneto máimo em cm
                
        if self.sist  == 'piso':                                                # Deslocamento máximo segundo: Specification for structural steel buildings (AISC 360-16)
            Ulim      = self.l / 350
        else:
            Ulim      = self.l / 250
                
        if Umax > Ulim:                                                         # Verificação do estado limite de serviço
            ver_ELS_AISC  = 'flecha excedida'
        else:
            ver_ELS_AISC  = 'OK'
            
#----------------------------------------------------------------------------#
# 3.7 - Limites de aplicabilidade do método AISC (2016)                      #
#----------------------------------------------------------------------------#             
       
        Lim1         = 'true'
        Lim2         = 'true'
        Lim3         = 'true'
        Lim_AISC     = 'true'
    
        if B_cond   == 'Não aplica':
            Lim1     = 'false'
            
        if self.bi == 0:
            r1       = round((self.p / self.d0),2)
            r2       = round((self.dg / self.d0),2)
            
            if r1 < 1.08 or r1 > 1.5:
                Lim2 = 'false'
            if r2 < 1.25 or r2 > 1.75:
                Lim3 = 'false'
        
        if Lim1 == 'false' or Lim2 == 'false' or Lim3 == 'false':
            Lim_AISC = 'false'

#----------------------------------------------------------------------------#
# 3.8 - Atribuições do método AISC ao objeto viga alveolar                   #
#----------------------------------------------------------------------------#                
        
        self.B_cond     = B_cond
        self.vsd_gross  = vsd_gross
        self.vsd_net    = vsd_net
        self.comp       = comp
        self.fle_AISC   = Umax
        self.flelim_AISC= Ulim
        self.Lim_AISC   = Lim_AISC
        self.Mebm       = Mebm
        self.Uv         = Uv
        self.Umax       = Umax
        self.B_pos_FLT  = B_pos_FLT
      
        return AISC_RRS,AISC_FMVa,AISC_FMVm,AISC_FMV,AISC_FMA,AISC_FMAC,    \
               AISC_FMACnet,AISC_FLT,ver_ELS_AISC
    
#----------------------------------------------------------------------------#
# 4.0 - Verificação por SCI No 100 and BS595                                 #
#----------------------------------------------------------------------------#     

    def SCI(self):  

# Valores de referência
    
        C_ver_FMPS   = 'OK'
        C_ver_RRS    = 'OK'
        C_ver_VCV    = 'OK'
        C_pos_FMPS   = -1
        C_pos_RRS    = -1
        C_pos_VCV    = -1
        C_pos_FMV    = -1
        C_pos_FLT    = -1

#----------------------------------------------------------------------------#
# 4.1 - Limites de aplicabilidade do método SCI No 100 and BS5950            #
#----------------------------------------------------------------------------#

        Lim1         = 'true'
        Lim2         = 'true'
        Lim_SCI      = 'true'

        if self.bi == 0:
            r1       = round((self.p / self.d0),2)
            r2       = round((self.dg / self.d0),2)
            
            if r1 < 1.08 or r1 > 1.5:
                Lim1 = 'false'
            if r2 < 1.25 or r2 > 1.75:
                Lim2 = 'false'
        
        if Lim1 == 'false' or Lim2 == 'false' or self.bi != 0:
            Lim_SCI  = 'false'

#----------------------------------------------------------------------------#
# 4.2 - Formação de Mecanismo Plástico Simples (FMPS)                        #
#----------------------------------------------------------------------------# 
    
        Mp_SCI        = self.At * self.y0 * 2 * self.fym                        # Momento de plastificação das seções T's   
        
        Mu            = np.max(self.vsd_ab)                                     # Momento atuante máximo na viga
        ind_FMPS      = np.argmax(self.vsd_ab)                                     # Indice do vetor que indica de maior esforço para falha por FMPS
        C_pos_FMPS    = self.sections[(ind_FMPS*4)+2]                           # Seção que ocorre falha por FMPS           
        Mpd           = Mp_SCI / self.coefseg                                   # Momento resistênte plástico das seções T's           
        
        if Mu > Mpd:                                                            # Verificação de falha por FMPS
            C_ver_FMPS = 'falha'
        
#----------------------------------------------------------------------------#
# 4.3 - Ruptura da região da solda (RRS)                                     #
#----------------------------------------------------------------------------#    

        vsdmax       = np.max(self.vsd)                                         # Cortante máximo atuante na viga
        VH           = self.vsd_mon * (self.p / (self.dg - 2*self.yl))          # Esforço cortante horizontal
        VHmax        = vsdmax * (self.p / (self.dg - 2*self.yb))                # Esforço cortante horizontal máximo
        Pvh          = 0.6 * self.fya * (0.9/self.coefseg) * self.bw * self.tw  # Resistência da região da solda entre aberturas
        Pvh_ap       = 0.6 * self.fya * (0.9/self.coefseg) * self.bwe * self.tw # Resistência da região da solda entre o apoio e a primeira abertura
        
        ind_vsdmax   = np.argmax(self.vsd)                                      # Indice correspondente ao cortante máximo do vetor de cortantes
        pos_vsdmax   = self.sections[ind_vsdmax]                                # Seção correspondente ao cortante máximo
       
        if VHmax > Pvh_ap and (pos_vsdmax == 0 or pos_vsdmax == self.l):
            C_ver_RRS= 'falha'
            C_pos_RRS= self.sections[ind_vsdmax]
        else:
            for x,i in enumerate(VH[:]):
                if i > Pvh:
                    C_ver_RRS= 'falha' 
                    C_pos_RRS= self.sections[(x+1)*4]         
            
#----------------------------------------------------------------------------#
# 4.4 - Verificação do Cisalhamento Vertical (VCV)                           #
#----------------------------------------------------------------------------# 
            
        Pvy           = 0.6*self.fya*(0.9/self.coefseg)*(2*self.ht*self.tw)     # Resistência da região da solda entre aberturas
        Pvy_ap        = (0.6 * self.fya * (0.9/self.coefseg) * 
                        ((self.dg-2*self.tf)*self.tw))                          # Resistência da região da solda na região de alma cheia

        ind_vsdmax    = np.argmax(self.vsd)                                     # Indice correspondente ao cortante máximo do vetor de cortantes
        pos_vsdmax    = self.sections[ind_vsdmax]                               # Seção correspondente ao cortante máximo       

        if pos_vsdmax == 0 or pos_vsdmax == self.l:            
            if vsdmax > Pvy_ap: 
                C_ver_VCV= 'falha'
                C_pos_VCV= self.sections[ind_vsdmax]
        else:
            for x,i in enumerate(self.vsd_mon[:]):
                if i > Pvy_ap:
                    C_ver_VCV= 'falha' 
                    C_pos_VCV= self.sections[(x+1)*4]
            for x,i in enumerate(self.vsd_ab[:]):
                if i > Pvy:
                    C_ver_VCV= 'falha' 
                    C_pos_VCV= self.sections[(x*4)+2]  
    
#----------------------------------------------------------------------------#
# 4.5 - Flambagem do Montante da Alma por Cisalhamento (FMA)                 #
#----------------------------------------------------------------------------#

        C1         = (5.097 + 0.1464 * (self.d0/self.tw) - 0.001740 * 
                     (self.d0/self.tw)**2)                                      # Constante de calculo 1 para vigas celulares
        C2         = (1.441 + 0.0625 * (self.d0/self.tw) - 0.000683 * 
                     (self.d0/self.tw)**2)                                      # Constante de calculo 1 para vigas celulares
        C3         = (3.645 + 0.0853 * (self.d0/self.tw) - 0.001080 * 
                     (self.d0/self.tw)**2)                                      # Constante de calculo 1 para vigas celulares
        
        Me_wpc     = (self.tw*(self.p - self.d0 + 0.564*self.d0)**2 * 
                      self.fya / 6)                                             # Web Post Capacity 
        
        Mmax         = Me_wpc * (C1 * (self.p/self.d0) - C2 * 
                       (self.p/self.d0)**2 - C3)                                # Momento máximo resistivo                  
        
        mdiff   = abs(self.msd_ab[1:] - self.msd_ab[:-1])                       # Momento atuante exatamente em cada poste de alma
        
        Vrh     = mdiff / self.deffec                                           # Cortante horizontal atuante no meio do poste da alma
        
        Mrh     = 0.9 * Vrh * (self.d0/2)                                       # Momento solicitante à 0.9 da altura do raio do alveolo
        
        C_ver_FMA     = 'não falha'                                             # Valor de refereência. será substituido caso haja falha
  
        for x,i in enumerate(Mrh[:]):
            if i > Mmax:
                C_ver_FMA     = 'falha'
                C_pos_FMA     = self.sections[(x+1)*4]
                self.C_pos_FMA= C_pos_FMA

#----------------------------------------------------------------------------#
# 4.6 - Formação do Mecanismo de Vierendeel (FMV)                            #
#----------------------------------------------------------------------------#

        Ti          = self.msd_ab / self.deffec                                 # Esforço horizontal na região das aberturas            
                                        
    # Esforços ba seção em teta = 25.84    
        
        teta_FMV    = 25                                                       

    # Capacidade segundo aproximação do método de Sahmel's

        tf_l        = self.tf / np.cos(np.radians(teta_FMV))
        l_l         = (((self.dg/2)-self.tf) / np.cos(np.radians(teta_FMV)) -\
                       (self.d0/2))
        At_l        = tf_l * self.bf + l_l * self.tw
        yl_l        = (tf_l*self.bf*(tf_l/2) + l_l*self.tw*((l_l/2)+tf_l))/At_l 
        ht_l        = tf_l + l_l
        
        Pu          = At_l * self.fym                                           # Capacidade axial da seção crítica na abertura     
        
        if At_l/(2*self.bf) < self.tf:
            Zx_l= (((self.tw*ht_l**2)/2) + ((self.bf*self.tf**2)/4) - \
                    ((ht_l*self.tf *self.tw)/2) - (((ht_l-self.tf)**2)/ \
                                                   (4*self.bf)))
        else:
            Zx_l= ((self.tw*(ht_l-self.tf)**2/4) + (self.bf*ht_l*self.tf/2) - \
                  ((self.bf**2*self.tf**2)/(4*self.tw)))    
        
        Mp_l        = Zx_l * self.fym                                           # Capacidade à flambagem da seção crítica na abertura

        Po          = np.cos(np.radians(teta_FMV)) * Ti - ((self.vsd_ab/2) * 
                      np.sin(np.radians(teta_FMV)))                             # Cortante no ponto 
        Mo          = (Ti * (yl_l - self.yl) + (self.vsd_ab/2) * 
                      ((self.dg/2) - yl_l)) * np.tan(np.radians(teta_FMV))      # Momento no ponto em questão

        R_FMV       = (Po / Pu) + (Mo / Mp_l)                                   # Relação de verificação
        
        C_ver_FMV   = 'não falha'                                               # Assume valor de não falha
  
        for x,i in enumerate(R_FMV[:]):                                         # Lopping de verificação de falha
            if i > (1/self.coefseg):
                C_ver_FMV       = 'falha'
                C_pos_FMV       = self.sections[(x+1)*4]
                self.C_pos_FMV  = C_pos_FMV 

#----------------------------------------------------------------------------#
# 4.7 - Flambagem Lateral com Toção (FLT)                                    #
#----------------------------------------------------------------------------#

# Verificação de tipo de seção (Compacta / semi-compacta / plástica / esbelta)

        epslon        = np.sqrt(275/self.fya)
        Fc            = np.max(Ti)
        r1            = Fc / ((self.dg-2*self.tf) * self.tw * self.fya)
        r2            = Fc / (self.Anet * self.fya)
        
        r_v           = ((self.dg-2*self.tf) - 3 * self.tw) / self.tw
        
        cl1           = 64 * epslon / (1 + 0.6 * r1)
        cl2           = 80 * epslon / (1 + r1)
        cl3           = 120 * epslon / (1 + 2 * r2) 
        
        if r_v <= cl1:
            tipo_sec  = 'compacta'
        elif r_v > cl1 and r_v <= cl2:
            tipo_sec  = 'semi-compacta'
        elif r_v > cl2 and r_v <= cl3:
            tipo_sec  = 'plastica'
        elif r_v > cl3:
            tipo_sec  = 'esbelta'
            
# Cálculo do fator equivalene unifome do momento resistente

        mfmax        = np.max(self.mf)                                          # Momento máximo atuante na viga
        tam_sec      = len(self.sections)                                       # Tamanho do vetor das seções analisadas pelas vigas alveolares
        ind_m2       = round(1/4 * tam_sec) - 1                                 # Posição no vetor de seções correspondente à 1/4 da viga
        m2           = self.mf[ind_m2]                                          # Momento atuante na posição 1/4 do comprimento da viga
        ind_m3       = round(2/4 * tam_sec)                                     # Posição no vetor das seções correspondente ao meio da viga 
        m3           = self.mf[ind_m3]                                          # Momento atuante no meio da viga
        ind_m4       = round(3/4 * tam_sec)                                     # Posição no vetor das seções corresponte à 3/4 da viga 
        m4           = self.mf[ind_m4]                                          # Momento atuante à 3/4 da viga
    
        mLT          = 0.2 + (0.15 * m2 + 0.5 * m3 + 0.15 * m4) / mfmax         # Fator equivalente para o momento uniforme
        
# Cálculo das esbeltezes

        lamb1        = self.lb / self.ry                                        # Comprimento de esbeltez
        lamb_L0      = 0.4 *(np.pi**2 * self.melast / self.fya)**0.5            # Comprimento equivalente de esbeltez
        
        if tipo_sec == 'compacta' or tipo_sec == 'semi-compacta':               # Razão 
            B_w      = 1
        elif tipo_sec == 'plastica':
            B_w      = self.Zx_vz / self.Wx_vz
        elif tipo_sec == 'esbelta':
            B_w      = 0.9 * self.Zx_vz / self.Wx_vz
        
        gama         = (1 - self.Iy_vz / self.Ix_vz)
        u            = (4 * self.Wx_vz**2 * gama / self.Anet**2 / 
                        self.deffec**2)**0.25
        x            = 0.566 * self.deffec * (self.Anet/self.J_vz)**0.5
        v            = 1 / (1 + 0.05 * (lamb1/x)**2)**0.25
        
        lamb_LT      = u * v * lamb1 * np.sqrt(B_w)                             # Comprimento equivalente de esbeltez
        
# Calculo da tensão escoamento equivalente

        alfa_LT      = 0.7                                                      # constante de Robertson
        ni_LT        = alfa_LT * (lamb_LT - lamb_L0) / 1000
        PE           = (np.pi**2 * self.melast) / lamb_LT**2
        phi_LT       = (self.fya + (ni_LT + 1) * PE) / 2 
 
        if lamb_LT < lamb_L0:                                                   # Definição da resistência à flambagem
            PB       = self.fya         
        else:    
            PB       = self.fya * PE / (phi_LT + (phi_LT**2 - 
                       PE*self.fya))**0.5 
        
# Calculo do momento resistente 

        if tipo_sec == 'compacta' or tipo_sec == 'semi-compacta':               # Definição do momento resistente 
            MB       = PB * self.Wx_vz
        elif tipo_sec == 'plastica':
            MB       = PB * self.Zx_vz
        elif tipo_sec == 'esbelta':
            MB       = 0.9 * PB * self.Zx_vz

# Verificação final

        C_ver_FLT     = 'não falha'     
        
        if mfmax > MB/mLT and self.lb != 0:                                                      # Verificação de falha
            C_ver_FLT     = 'falha'
            ind_FLT       = np.argmax(self.mf)
            C_pos_FLT     = self.sections[(ind_FLT)]
            self.C_pos_FLT= C_pos_FLT
            
#----------------------------------------------------------------------------#
# 4.8 - Verificação do Deslocamento Máximo (ELS)                             #
#----------------------------------------------------------------------------#

# Deflexão segundo a SCI e com acrescimo de 25% segundo Pachor et al. (2014)

        if self.tipo == 'distribuida':
            w        = self.q/100
            U        = (5/384) * w * self.l**4 / self.melast / self.Ix_vz
        elif self.tipo == 'concentrada': 
            U        = (self.q * self.l**3 / (48 * self.melast * self.Ix_vz))
                    
        Umax          = 1.25 * U
                
        if self.sist  == 'piso':
            Ulim      = self.l / 350
        else:
            Ulim      = self.l / 250
                
        if Umax > Ulim:
            C_ver_ELS  = 'flecha excedida'
        else:
            C_ver_ELS  = 'OK'

#----------------------------------------------------------------------------#
# 4.9 - Atribuições do método SCI ao objeto viga alveolar                    #
#----------------------------------------------------------------------------#
        
        self.C_pos_FMPS     = C_pos_FMPS
        self.C_pos_RRS      = C_pos_RRS
        self.C_pos_VCV      = C_pos_VCV
        self.C_pos_FLT      = C_pos_FLT
        self.Lim_SCI        = Lim_SCI
        self.C_ver_ELS      = C_ver_ELS
        self.SCI_Umax       = Umax
        self.SCI_Ulim       = Ulim

        return C_ver_FMPS,C_ver_RRS,C_ver_VCV,C_ver_FMA,C_ver_FMV,C_ver_FLT,\
               C_ver_ELS
               
#----------------------------------------------------------------------------#
# 5.0 - Verificação por Grilo et al. (2018)                                  #
#----------------------------------------------------------------------------#     

    def grilo(self):  
        
        D_ver_FMA   = 'OK'
        D_pos_FMA   = -1
        
        # Extraido de Veríssimo et al. (2012)
        
        D_ver_FMPS = 0
        D_ver_RRS  = 0
        D_ver_FLT  = 0 
        pos_FMPS   = -1                                                               
        pos_RRS    = -1                                                                                                                             
        pos_FLT    = -1

#----------------------------------------------------------------------------#
# 5.1 - Limites de aplicabilidade do método GRILO et al. (2018)              #
#----------------------------------------------------------------------------#

        Lim1         = 'true'
        Lim2         = 'true'
        Lim3         = 'true'
        Lim_grilo    = 'true'

        if self.bi == 0:
            r1       = round((self.p / self.d0),1)
            r2       = round((self.d0 / self.dg),1)
            LM       = 0.5 * np.sqrt(self.p**2 - self.d0**2)
            r3       = LM * np.sqrt(12) / self.tw
            
            if r1 < 1.1 or r1 > 1.5:
                Lim1 = 'false'
            if r2 < 0.55 or r2 > 0.8:
                Lim2 = 'false'
            if r3 < 10 or r3 > 200:
                Lim3 = 'false'
        
        if Lim1 == 'false' or Lim2 == 'false' or Lim3 == 'false' or \
           self.bi != 0:
            Lim_grilo= 'false'
            
#----------------------------------------------------------------------------#
# 5.2 - Flambagem do Montante da Alma por Cisalhamento (FMA)                 #
#----------------------------------------------------------------------------#            

# Calculo do beta

        if self.p / self.d0 < 1.2:
            beta     = 1.198 - 0.42 * (self.d0/self.dg) + self.p / (5*self.d0)
        elif self.p / self.d0 >= 1.2: 
            beta     = 1.838 - 0.42 * (self.d0/self.dg) - self.p / (3*self.d0)
            
# Altura da plastificação 

        yp           = (self.d0 / 2) * (0.445 * (self.p/self.d0)**3 - 2.578 * 
                       (self.p / self.d0)**2 + 4.770 * (self.p/self.d0)-2.475)
        
# Largura do montante da alma a uma altura yp

        bm_yp        = self.p - self.d0 * np.sqrt(1-(4*yp**2/self.d0**2))
        
# Esforço cortante resistente de plastificação em uma seção à altura yp

        Vhp          = (beta * self.fya * self.tw * bm_yp**2 / 
                        np.sqrt(3*bm_yp**2 + 16*yp**2))
        # Vhp_ap       = (beta * self.fya * self.tw * self.bwe**2 / 
        #                 np.sqrt(3*self.bwe**2 + 16*yp**2))
        
# Relações de cálculo para abaco

        k1           = round(self.d0/self.dg,1)        
        k2           = round(self.p/self.d0,1)       
        
# Valores do abaco 
       
        ax, ay       = np.meshgrid(self.op_a.columns,self.op_a.index)
        at           = interp2d(ax,ay,self.op_a.values)
        a            = at(k1,k2)
        
        bx, by       = np.meshgrid(self.op_b.columns,self.op_b.index)
        bt           = interp2d(bx,by,self.op_b.values)
        b            = bt(k1,k2)
        
        cx, cy       = np.meshgrid(self.op_c.columns,self.op_c.index)
        ct           = interp2d(cx,cy,self.op_c.values)
        c            = ct(k1,k2)
        
        dx, dy       = np.meshgrid(self.op_d.columns,self.op_d.index)
        dt           = interp2d(dx,dy,self.op_d.values)
        d            = dt(k1,k2)
        
        ex, ey       = np.meshgrid(self.op_e.columns,self.op_e.index)
        et           = interp2d(ex,ey,self.op_e.values)
        e            = et(k1,k2)
        
        self.c = c
        self.d = d
        
# Indice de esbeltez reduzida 

        lamb_ma_0    = np.sqrt((3*(self.p**2 - self.d0**2)*self.fya)/
                               (np.pi**2*self.tw**2*self.melast))
        
# Calculo do fator de redução X 

        if lamb_ma_0 >= 1.0:
            X        = a / lamb_ma_0**b
        else:  
            X        = c * d**(lamb_ma_0**e)
            
        if X > 1:
            X = 1          

# Verificação de FMV

        if c <= 0.4 or d <= 0.4:
            D_ver_FMA ='Indicação de falha por FMV'
        else: D_ver_FMA ='OK'
        
# Calculo do esfoço horizontal solicitante 
         
        Vhs             = self.vsd_mon * (self.p/(2*self.y0))

# Esforço cortante resistente à FMA

        Vrk             = X * Vhp

# Verificação de falha    

        if D_ver_FMA == 'OK':
            for x,i in enumerate(Vhs[:]):
                if i > Vrk:
                    D_ver_FMA = 'Falha por FMA'
                    D_pos_FMA = self.sections[(x+1)*4]
                    self.D_pos_FMA  = D_pos_FMA

#----------------------------------------------------------------------------#
# 5.3 - Verificação de Flecha Máxima (ELS)                                   #
#----------------------------------------------------------------------------#
        
        if self.tipo == 'distribuida':
            w         = self.q/100
            Um        = (5/384) * w* self.l**4 / self.melast / self.Ie
            Uad       = w * self.l**2 / 8 / self.G / self.Ae
            Uv        = 0
        elif self.tipo == 'concentrada': 
            alfa      = (self.Ae/(self.tw*self.Ie)) * ((self.bf*self.dg**2/8) 
                        - ((self.dg-2*self.tf)**2/8) * (self.bf-self.tw))
            Uad       = (alfa / self.Ae / self.G) * (self.q * self.l / 4)
            Uv        = (self.q * self.l**3 / (48 * self.melast * self.Ie))
            Um        = 0
                    
        Umax          = (Um + Uv + Uad)
                
        if self.sist  == 'piso':
            Ulim      = self.l / 350
        else:
            Ulim      = self.l / 250
                
        if Umax > Ulim:
            D_ver_ELS = 'flecha excedida'
        else:
            D_ver_ELS = 'OK'

# AS DEMAIS VERIFICAÇÕES SÃO DO MÉTODO DO VERÍSSIMO et al.(2012)

#----------------------------------------------------------------------------#
# 5.4 - Formação de Mecanismo Plástico Simples (FMPS)                        #
#----------------------------------------------------------------------------# 
     
        c             = self.y0*self.ya*self.bw*self.At/(2*self.Ix_t)           # Constante que multiplica o esforço cortante (DELESQUES, 1969)
        distesf       = self.vsd*c + self.mf                                    # Distância entre a diferença dos esforços cortante equivalente e momento fletor
       
        distaux       = np.zeros(int(self.n/2)+2)                                        
        distaux[0]    = distesf[2]                                                   
        distaux[1]    = distesf[6]                                                   
        distaux[2:]   = distesf[10::4][:int(self.n/2)]
       
        esfsec        = np.max(distaux)                                         # Maior distância entre os esforços encontrados ao longo da viga
        ind_FMPS      = np.argmax(distaux)                                      # Indice do vetor que indica de maior esforço para falha por FMPS
        pos_FMPS      = self.sections[(ind_FMPS*4)+2]                           # Seção que ocorre falha por FMPS           
        Mpl           = self.fym*self.Zx_vz                                     # Momento de plastificação das seções T's
        Mrd           = Mpl * self.coefseg                                      # Momento resistente (Mpl x Yal) por norma da seção T
        Mesl          = esfsec                                                  # Momento atuante máximo para falha por FMPS
        
        if Mesl > Mrd:                                                          # Condição de falha à FMPS
            D_ver_FMPS= 'falha'

#----------------------------------------------------------------------------#
# 5.5 - Ruptura na Regição da Solda (RRS)                                    #
#----------------------------------------------------------------------------# 
        
        vsdmax        = np.max(self.vsd_mon)                                    # Cortante máximo atuante nos montante de alma na viga
        
        vrk1          = (4/(3*np.sqrt(3))) * (self.bw*self.tw*self.y0*
                         self.fya/self.p)                                       # Esforço cortante resistente característico ao E.M.A.F.
        
        Vrd1          = vrk1 / self.coefseg                                     # Cortante de cálculo ao E.M.A.F.
        
        ind_vsdmax    = np.argmax(vsdmax)                                       # Indice correspondente ao cortante máximo do vetor de cortantes
        pos_vsdmax    = self.sections[(ind_vsdmax+1)*4]                         # Seção correspondente ao cortante máximo
        
        if vsdmax > Vrd1:
            D_ver_RRS= 'falha'
            pos_RRS   = pos_vsdmax
                         
#----------------------------------------------------------------------------#
# 5.6 - Flambagem Lateral com Torção (FLT)                                   #
#----------------------------------------------------------------------------#        
     
        mfmax        = np.max(self.mf)                                          # Momento máximo atuante na viga
        tam_sec      = len(self.sections)                                       # Tamanho do vetor das seções analisadas pelas vigas alveolares
        ind_ma       = round(1/4 * tam_sec) - 1                                 # Posição no vetor de seções correspondente à 1/4 da viga
        ma           = self.mf[ind_ma]                                          # Momento atuante na posição 1/4 do comprimento da viga
        ind_mb       = round(2/4 * tam_sec)                                     # Posição no vetor das seções correspondente ao meio da viga 
        mb           = self.mf[ind_mb]                                          # Momento atuante no meio da viga
        ind_mc       = round(3/4 * tam_sec)                                     # Posição no vetor das seções corresponte à 3/4 da viga 
        mc           = self.mf[ind_mc]                                          # Momento atuante à 3/4 da viga
        cb           = 12.5*mfmax / (2.5*mfmax + 3*ma + 4*mb + 3*mc)            # Constante de cálculo (Usado para determinar o momento característico)
        beta1        = 0.7*self.fym*self.Wx_vz/(self.melast*self.J_vz)          # Constante de cálculo
        Lp           = (1.76*self.ry_vz*np.sqrt(self.melast/self.fym))          # Comprimento de plastificação                                                                
        Lrcor        = (1.66*np.sqrt(self.Iy_vz*self.J_vz)/
                       (self.J_vz*beta1))* (np.sqrt(1+np.sqrt(1+
                       (27*self.Cw*beta1**2/self.Iy_vz))))                      # Comprimento de escoamento
        Mrcor        = (0.31*self.melast/Lrcor**2) * (np.sqrt(self.Iy_vz*
                       (1000*self.Cw + 39*self.J_vz*self.lb**2)))               # Momento de início de escoamento
        
        if self.lb > Lrcor :                                                    # Determinação do momento resistente caracteristico   
            Mrk      = ((cb*np.pi**2*self.melast*self.Iy_vz/self.lb**2) *
                       (np.sqrt((self.Cw/self.Iy_vz) * 
                       (1+0.39*self.J_vz*self.lb**2/self.Cw))))
        elif Lp < self.lb and self.lb <= Lrcor:
            Mrk      = (cb * (0.9 * Mpl - (0.9 * Mpl - Mrcor) * 
                       ((self.lb - Lp)/(Lrcor - Lp))))
        else:
            Mrk      = 0.9 * Mpl
        
        ind_FLT      = np.argmax(self.mf)                                       # Indice correspondente ao momento máximo do vetor de cortantes  
        pos_FLT      = self.sections[ind_FLT]                                   # Posição que ocorre o momento máximo
        
        Msd_FLT      = Mrk / self.coefseg                                       # Momento resistente à FLT
        
        if mfmax > Msd_FLT:                                                     # Verificação de falha
            D_ver_FLT  = 'falha'
        
#----------------------------------------------------------------------------#
# 5.4 - Atribuições do método Grilo et al. (2018) ao objeto viga alveolar    #
#----------------------------------------------------------------------------#

        self.D_ver_ELS    = D_ver_ELS
        self.fle_grilo    = Umax
        self.flelim_grilo = Ulim 
        self.D_pos_FMA    = D_pos_FMA
        self.pos_FMPS     = pos_FMPS
        self.pos_RRS      = pos_RRS
        self.pos_FLT      = pos_FLT
        
        return Lim_grilo,D_ver_FMA,D_ver_ELS,D_ver_FMPS,D_ver_RRS,D_ver_FLT

#----------------------------------------------------------------------------#
# 6.0 - Anexo N / Eurocode 3 (ENV 1993-1-1:1992/A2:1998)                     #
#----------------------------------------------------------------------------#

    def annexN(self): 

# Valores de referência
        
        E_ver_FMPS   = 'OK'
        E_ver_RRS    = 'OK'
        E_ver_VCV    = 'OK'
        E_pos_FMPS   = -1
        E_pos_RRS    = -1
        E_pos_VCV    = -1
        E_pos_FMV    = -1
        
#----------------------------------------------------------------------------#
# 6.1 - Limites de aplicabilidade do método                                  #
#----------------------------------------------------------------------------#        
        
        if self.bi == 0 and (self.bw < 0.25*self.d0 or self.bw > 0.55*self.d0):  # Apenas para vigas celulare
            Lim_annexN= 'false'
        elif self.bi != 0:
            Lim_annexN= 'false'
        else:
            Lim_annexN   = 'true'

#----------------------------------------------------------------------------#
# 6.2 - Formação de Mecanismo Plástico Simples (FMPS)                        #
#----------------------------------------------------------------------------#
    
        Vplrd         = self.Ag_vz * self.fym / np.sqrt(3) / self.coefseg       # Resistência de projeto plástica ao corte  
        vsdmax        = np.max(self.vsd)                                        # Máximo cortante atuante na viga
        mfmax         = np.max(self.mf)                                         # Máximo momento atuante na viga
        RO            = ((2 * vsdmax / Vplrd) - 1)**2                           # Constante de cálculo  
        ind_FMPS      = np.argmax(self.mf)                                      # Indice do vetor que indica de maior esforço para falha por FMPS
        E_pos_FMPS    = self.sections[ind_FMPS]                                 # Seção que ocorre falha por FMPS 
         
        if vsdmax > 0.5 * Vplrd:
            Mvrd      = ((self.Zx_vz - (RO * self.Anet**2 / 4 / self.tw)) *  
                             self.fym) / self.coefseg                           # Momento resistente
        else:
            Mvrd      = (self.Zx_vz * self.fym) / self.coefseg                  # Momento resistente
            
        if mfmax > Mvrd:                                                        # Verificação de falha por FMPS  
            E_ver_FMPS= 'falha' 

#----------------------------------------------------------------------------#
# 6.3 - Formação de Mecanismo de Vierendeel (FMV)                            #
#----------------------------------------------------------------------------#

        Ti          = self.msd_ab / self.deffec                                 # Esforço horizontal na região das aberturas            
                                        
    # Esforços ba seção em teta = 25.84    
        
        teta_FMV    = 25                                                       

    # Capacidade segundo aproximação do método de Sahmel's

        tf_l        = self.tf / np.cos(np.radians(teta_FMV))
        l_l         = (((self.dg/2)-self.tf) / np.cos(np.radians(teta_FMV)) -\
                       (self.d0/2))
        At_l        = tf_l * self.bf + l_l * self.tw
        yl_l        = (tf_l*self.bf*(tf_l/2) + l_l*self.tw*((l_l/2)+tf_l))/At_l 
        ht_l        = tf_l + l_l
        
        Pu          = At_l * self.fym                                           # Capacidade axial da seção crítica na abertura     
        
        if At_l/(2*self.bf) < self.tf:
            Zx_l= (((self.tw*ht_l**2)/2) + ((self.bf*self.tf**2)/4) - \
                    ((ht_l*self.tf *self.tw)/2) - (((ht_l-self.tf)**2)/ \
                                                   (4*self.bf)))
        else:
            Zx_l= ((self.tw*(ht_l-self.tf)**2/4) + (self.bf*ht_l*self.tf/2) - \
                  ((self.bf**2*self.tf**2)/(4*self.tw)))    
        
        Mp_l        = Zx_l * self.fya                                           # Capacidade à flambagem da seção crítica na abertura

        Po          = np.cos(np.radians(teta_FMV)) * Ti - ((self.vsd_ab/2) * 
                      np.sin(np.radians(teta_FMV)))                             # Cortante no ponto 
        Mo          = (Ti * (yl_l - self.yl) + (self.vsd_ab/2) * 
                      ((self.dg/2) - yl_l)) * np.tan(np.radians(teta_FMV))      # Momento no ponto em questão

        R_FMV       = (Po / Pu) + (Mo / Mp_l)                                   # Relação de verificação
        
        E_ver_FMV       = 'não falha'
        
        for x,i in enumerate(R_FMV[:]):                                         # Lopping de verificação de falha
            if i > 1:
                E_ver_FMV       = 'falha'
                E_pos_FMV       = self.sections[(x+1)*4]
                self.E_pos_FMV  = E_pos_FMV  

#----------------------------------------------------------------------------#
# 6.4 - Flambagem do Montante da Alma por Cisalhamento (FMA)                 #
#----------------------------------------------------------------------------#

        C1         = (5.097 + 0.1464 * (self.d0/self.tw) - 0.001740 * 
                     (self.d0/self.tw)**2)                                      # Constante de calculo 1 para vigas celulares
        C2         = (1.441 + 0.0625 * (self.d0/self.tw) - 0.000683 * 
                     (self.d0/self.tw)**2)                                      # Constante de calculo 1 para vigas celulares
        C3         = (3.645 + 0.0853 * (self.d0/self.tw) - 0.001080 * 
                     (self.d0/self.tw)**2)                                      # Constante de calculo 1 para vigas celulares
        
        Me_wpc     = (self.tw*(self.p - self.d0 + 0.564*self.d0)**2 * 
                      self.fya / 6)                                             # Web Post Capacity 
        
        Mmax         = Me_wpc * (C1 * (self.p/self.d0) - C2 * 
                       (self.p/self.d0)**2 - C3)                                # Momento máximo resistivo                  
        
        mdiff   = abs(self.msd_ab[1:] - self.msd_ab[:-1])                       # Momento atuante exatamente em cada poste de alma
        
        Vrh     = mdiff / self.deffec                                           # Cortante horizontal atuante no meio do poste da alma
        
        Mrh     = 0.9 * Vrh * (self.d0/2)                                       # Momento solicitante à 0.9 da altura do raio do alveolo
       
        E_ver_FMA  = 'não falha'                                                # Valor de refereência. será substituido caso haja falha  
 
        for x,i in enumerate(Mrh[:]):                                           # Verificação de falha ao longo da viga
             if i > Mmax:
                 E_ver_FMA     = 'falha'
                 E_pos_FMA     = self.sections[(x+1)*4]
                 self.E_pos_FMA= E_pos_FMA

#----------------------------------------------------------------------------#
# 6.5 - Verificação do Cisalhamento Vertical (VCV)                           #
#----------------------------------------------------------------------------# 
 
        Pvy           = self.fya * (1/self.coefseg) * self.Anet / np.sqrt(3)    # Resistência da região da solda entre aberturas
        Pvy_ap        = self.fya * (1/self.coefseg) * self.Ag / np.sqrt(3)      # Resistência da região da solda entre o apoio e a primeira abertura

        ind_vsdmax    = np.argmax(self.vsd)                                     # Indice correspondente ao cortante máximo do vetor de cortantes
        pos_vsdmax    = self.sections[ind_vsdmax]                               # Seção correspondente ao cortante máximo       

        if pos_vsdmax == 0 or pos_vsdmax == self.l:            
            if vsdmax > Pvy_ap: 
                E_ver_VCV= 'falha'
                E_pos_VCV= self.sections[ind_vsdmax]
        else:
            for x,i in enumerate(self.vsd_mon[:]):
                if i > Pvy_ap:
                    E_ver_VCV= 'falha' 
                    E_pos_VCV= self.sections[(x+1)*4]
            for x,i in enumerate(self.vsd_ab[:]):
                if i > Pvy:
                    E_ver_VCV= 'falha' 
                    E_pos_VCV= self.sections[(x*4)+2]

#----------------------------------------------------------------------------#
# 6.6 - Ruptura na região da solda (RRS)                                     #
#----------------------------------------------------------------------------#

        vsdmax        = np.max(self.vsd)                                        # Cortante máximo atuante na viga
        VH            = self.vsd_mon * (self.p / (self.dg - 2*self.yb))         # Esforço cortante horizontal
        VHmax         = vsdmax * (self.p / (self.dg - 2*self.yb))               # Esforço cortante horizontal máximo
        Pvh           = self.fya*self.bw*self.tw/(np.sqrt(3)*self.coefseg)      # Resistência da região da solda entre aberturas
        Pvh_ap        = self.fya*self.bwe*self.tw/(np.sqrt(3)*self.coefseg)     # Resistência da região da solda entre o apoio e a primeira abertura
       
        ind_vsdmax    = np.argmax(self.vsd)                                     # Indice correspondente ao cortante máximo do vetor de cortantes
        pos_vsdmax    = self.sections[ind_vsdmax]                               # Seção correspondente ao cortante máximo
             
        if VHmax > Pvh_ap and (pos_vsdmax == 0 or pos_vsdmax == self.l):
            E_ver_RRS= 'falha'
            E_pos_RRS= self.sections[ind_vsdmax]
        else:
            for x,i in enumerate(VH[:]):
                if i > Pvh:
                    E_ver_RRS= 'falha' 
                    E_pos_RRS= self.sections[(x+1)*4] 

#----------------------------------------------------------------------------#
# 6.7 - Flambagem lateral com torção (FLT)                                   #
#----------------------------------------------------------------------------#        

#----------------------------------------------------------------------------#
# 6.8 - Verificação do Deslocamento Máximo (ELS)                             #
#----------------------------------------------------------------------------#

# Deflexão segundo a SCI e com acrescimo de 25% segundo Pachor et al. (2014)

        if self.tipo == 'distribuida':
            w         = self.q/100
            Um        = (5/384) * w * self.l**4 / self.melast / self.Ix_vz
        elif self.tipo == 'concentrada': 
            Um        = (self.q * self.l**3 / (48 * self.melast * self.Ix_vz))
                    
        Umax          = 1.25 * Um
                
        if self.sist  == 'piso':
            Ulim      = self.l / 350
        else:
            Ulim      = self.l / 250
                
        if Umax > Ulim:
            E_ver_ELS  = 'flecha excedida'
        else:
            E_ver_ELS  = 'OK'

#----------------------------------------------------------------------------#
# 6.9 - Atribuições do método SCI ao objeto viga alveolar                    #
#----------------------------------------------------------------------------#
        
        self.E_pos_FMPS      = E_pos_FMPS
        self.E_pos_RRS       = E_pos_RRS
        self.E_pos_VCV       = E_pos_VCV
        self.E_ver_ELS       = E_ver_ELS
        self.Lim_annexN      = Lim_annexN
        self.annexN_Umax     = Umax
        self.annexN_Ulim     = Ulim      
    
        return E_ver_FMPS,E_ver_FMV,E_ver_FMA,E_ver_VCV,E_ver_RRS,E_ver_ELS
               
#----------------------------------------------------------------------------#
# 7.0 - Viga Original segundo a NBR 8800:2008                                #
#----------------------------------------------------------------------------#

    def original(self):  
        
        O_ver_FLM  = 'ok'
        O_ver_FLA  = 'ok'
        O_ver_FCR  = 'OK'
        O_ver_FLT  = 'OK'
        O_pos_FLM  = -1
        O_pos_FLA  = -1
        O_pos_FCR  = -1
        O_pos_FLT  = -1 

# Momentos de cálculo

        My_a         = self.Wx_O * self.fya                                     # Momento de início de plastificação da alma
        My_m         = self.Wx_O * self.fym                                     # Momento de início de plastificação da mesa
        Mp_a         = self.Zx_O * self.fya                                     # Momento de plastificação da alma
        Mp_m         = self.Zx_O * self.fym                                     # Momento de plastificação da mesa
        
#----------------------------------------------------------------------------#
# 7.1 - Dimensionamento a flexão - Flambagem Local da Mesa (FLM)             #
#----------------------------------------------------------------------------#         

# Determinação de seção compacta, semi-compacta ou esbelta

        lamb_bm    = self.bf / (2 * self.tf)
        lamb_pm    = 0.38* np.sqrt(self.melast / self.fym)
        kc         = 4 / np.sqrt(self.d/self.tw)
        C          = 0.83 
        lamb_rm    = C * np.sqrt(self.melast * kc / (0.7 * self.fym))                
        
        if lamb_bm <= lamb_pm:
            tipo_m = 'compacta'
        elif lamb_bm > lamb_pm and lamb_bm <= lamb_rm:
            tipo_m = 'semi-compacta'
        elif lamb_bm > lamb_rm:
            tipo_m = 'esbelta'
        
# Momento resistente de projeto

        Mr_m  = self.Wx_O * self.fym * 0.7
        if tipo_m == 'compacta':  
            Mn_m  = Mr_m
        elif tipo_m == 'semi-compacta':
            Mr_m  = Mp_m - (((lamb_bm-lamb_pm)/(lamb_rm-lamb_pm)) * 
                            (Mp_m - Mr_m))
            Mn_m  = Mr_m 
        elif tipo_m == 'esbelta':
            Mr_m  = 0.69 * self.melast * self.Wx_O / lamb_bm**2 
            Mn_m  = Mr_m
        
        Mrk_m     = Mn_m / self.coefseg
        
        if Mrk_m > 1.5 * self.Wx_O * self.fym / self.coefseg:
            Mrk_m = 1.5 * self.Wx_O * self.fym / self.coefseg 
            
# Verificação de resistencia

        Mmax      = np.max(self.mf)
        
        O_ind_FLM      = np.argmax(self.mf)                                     # Indice correspondente ao momento máximo do vetor de cortantes  
        O_pos_FLM      = self.sections[O_ind_FLM]                               # Posição que ocorre o momento máximo
        
        if Mmax > Mrk_m:                                                        # Verificação de falha
            O_ver_FLM  = 'falha'
            
#----------------------------------------------------------------------------#
# 7.2 - Dimensionamento a flexão - Flambagem Local da Alma (FLA)             #
#----------------------------------------------------------------------------# 

# Determinação de seção compacta, semi-compacta ou esbelta

        lamb_ba    = (self.d-2*self.tf) / self.tw
        D          = 3.76
        lamb_pa    = D * np.sqrt(self.melast/self.fya) 
        lamb_ra    = 5.70 * np.sqrt(self.melast/self.fya) 
        
        if lamb_ba <= lamb_pa:
            tipo_a = 'compacta'
        elif lamb_ba > lamb_pa and lamb_ba <= lamb_ra:
            tipo_a = 'semi-compacta'
        elif lamb_ba > lamb_ra:
            tipo_a = 'esbelta'
            
# Momento resistente de projeto

        ar        = (self.d-2*self.tf) * self.tw / (self.bw*self.tf)
        k         = 1 - (ar/(1200+300*ar)) * (lamb_ba - lamb_ra)        
        Me1       = self.Wx_O * self.fya        
        Me2       = self.Wx_O * self.fya * k
        
        if tipo_a == 'compacta':  
            Mr_a  = self.Wx_O * self.fym * 0.7
            Mn_a  = Mr_a
        elif tipo_a == 'semi-compacta':
            Mr_a  = Mp_a - (((lamb_ba-lamb_pa)/(lamb_ra-lamb_pa)) * 
                            (Mp_a - Mr_a))
            Mn_a  = Mr_a 
        elif tipo_a == 'esbelta':
            Mr_a  = min(Me1,Me2)
            Mn_a  = Mr_a

        Mrk_a     = Mn_a / self.coefseg          
        
        if Mrk_a > 1.5 * self.Wx_O * self.fya / self.coefseg:
            Mrk_a = 1.5 * self.Wx_O * self.fya / self.coefseg 
            
# Verificação de resistencia
        
        O_ind_FLA      = np.argmax(self.mf)                                     # Indice correspondente ao momento máximo do vetor de cortantes  
        O_pos_FLA      = self.sections[O_ind_FLA]                               # Posição que ocorre o momento máximo
        
        if Mmax > Mrk_a:                                                        # Verificação de falha
            O_ver_FLA  = 'falha'            

#---------------------------------------------------------------------------#
# 7.3 - Dimensionamento a flexão - Flambagem ateral com torção (FLT)        #
#---------------------------------------------------------------------------# 
                                                                                                      
        tam_sec      = len(self.sections)                                       # Tamanho do vetor das seções analisadas pelas vigas alveolares
        ind_ma       = round(1/4 * tam_sec) - 1                                 # Posição no vetor de seções correspondente à 1/4 da viga
        ma           = self.mf[ind_ma]                                          # Momento atuante na posição 1/4 do comprimento da viga
        ind_mb       = round(2/4 * tam_sec)                                     # Posição no vetor das seções correspondente ao meio da viga 
        mb           = self.mf[ind_mb]                                          # Momento atuante no meio da viga
        ind_mc       = round(3/4 * tam_sec)                                     # Posição no vetor das seções corresponte à 3/4 da viga 
        mc           = self.mf[ind_mc]                                          # Momento atuante à 3/4 da viga
        cb           = 12.5*Mmax / (2.5*Mmax + 3*ma + 4*mb + 3*mc)    

        lbp          = 1.76 * self.ry_O * np.sqrt(self.melast/self.fya)
        beta1        = 0.7 * self.fya * self.Wy_O / (self.melast * self.J_O)
        lbr          = ((1.38 * np.sqrt(self.Iy_O*self.J_O) / 
                       (self.J_O*beta1)) * (np.sqrt(1 + np.sqrt(1 + 
                       (27*self.Cw_O*beta1**2/self.Iy_O)))))

        if self.lb <= lbp:
            Mn       = (self.Zx_O * self.fya / self.coefseg) / self.coefseg
        elif self.lb > lbp and self.lb <= lbr:
            Mn       = (cb * (Mp_a - (Mp_a-(0.7 * self.fya * self.Wx_O)) * 
                       ((self.lb - lbp)/(lbr - lbp)))) / self.coefseg
        elif self.lb > lbr:
            Mn       = ((cb * (np.pi**2 * self.melast * self.Iy_O / 
                        self.lb**2) * (np.sqrt((self.Cw_O/self.Iy_O) * 
                       (1 + 0.039 * self.J_O * self.lb**2 / self.Cw_O))))/ 
                        self.coefseg)
            
        for x,i in enumerate(self.mf[:]):
                if i > Mn:
                    O_ver_FLT  = 'falha'
                    O_pos_FLT  = self.sections[x]
                    self.O_pos_FLT  = O_pos_FLT 
        
#---------------------------------------------------------------------------#
# 7.4 - Dimensionamento ao cortante - Força cortante resistent (FCR)        #
#---------------------------------------------------------------------------# 

        kv              = 5                                                     # constante que leva em conta sem enrijecedor transversal
        lamb_c          = (self.d - 2*self.tf) / self.tw 
        lamb_cp         = 1.10 * np.sqrt(kv * self.melast / self.fya)
        lamb_cr         = 1.37 * np.sqrt(kv * self.melast / self.fya)
        
        Aw              = self.d * self.tw
        Vpl             = 0.6 * Aw * self.fya
        
        if lamb_c <= lamb_cp:
            Vrd    = Vpl/self.coefseg
        elif lamb_c > lamb_cp and lamb_c <= lamb_cr:
            Vrd    = (lamb_cp/lamb_c) * (Vpl/self.coefseg)
        elif lamb_c > lamb_cr:
            Vrd    = 1.24 * (lamb_cp/lamb_c)**2 * (Vpl/self.coefseg)
            
        for x,i in enumerate(self.vsd[:]):
                if i > Vrd:
                    O_ver_FCR  = 'falha'
                    O_pos_FCR = self.sections[x]
                    self.O_pos_FCR  = O_pos_FCR
        
#---------------------------------------------------------------------------#
# 7.5 - Verificação do estado limite de serviço (ELS) - Flexão limite       #
#---------------------------------------------------------------------------# 

        if self.tipo == 'distribuida':
            w         = self.q/100
            Um        = (5/384) * w* self.l**4 / self.melast / self.Ix_O
            Uad       = w * self.l**2 / 8 / self.G / self.Ag_O
            Uv        = 0
        elif self.tipo == 'concentrada': 
            alfa      = (self.Ag_O/(self.tw*self.Ix_O)) * ((self.bf*self.d**2/8) 
                        - ((self.d-2*self.tf)**2/8) * (self.bf-self.tw))
            Uad       = (alfa / self.Ag_O / self.G) * (self.q * self.l / 4)
            Uv        = (self.q * self.l**3 / (48 * self.melast * self.Ix_O))
            Um        = 0
                    
        Umax          = (Um + Uv + Uad)
                
        if self.sist  == 'piso':
            Ulim      = self.l / 350
        else:
            Ulim      = self.l / 250
                
        if Umax > Ulim:
            O_ver_ELS  = 'flecha excedida'
        else:
            O_ver_ELS  = 'OK' 

#---------------------------------------------------------------------------#
# 7.6 - Atribuições da verificação da viga de perfil de alma cheia pela NBR #
#---------------------------------------------------------------------------#            
        
        self.O_pos_FLM   = O_pos_FLM        
        self.O_pos_FLA   = O_pos_FLA
        self.Ori_fle     = Umax
        self.Ori_flelim  = Ulim
        self.Mrk_m       = Mrk_m
        self.Mmax        = Mmax

        return O_ver_FLM,O_ver_FLA,O_ver_FCR,O_ver_FLT,O_ver_ELS 
