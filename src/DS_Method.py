#!/bin/env python3

"""
Created on Mon Apr  16 22:58:30 2018

@author: pedro.moresco
v0.0.1
This algorithm is the adaptation for an API lybrary for calculate the stress
in structures by using the direct stiffness method.
Next tasks:
--> add a function to use txt archives
-->add a port to plot the structure
-->prepare API! - ready to go and tested Tue Apr 17 23:16
-->add the indication of which the bar is in compression or traction
-->translate to english
-->add docstrings and help
-->pass nodes as a list
"""

import numpy as np

def stiff_matrices(bars, nodes, area, mom_e):
    forcas_index = []
    fi = 0
    for ff in range(len(nodes)):
        forcas_index.append((fi, fi+1))
        fi += 2
    mrl = np.zeros((4, 4, len(bars)))
    # calcula o comprimento das barras
    comp_barras = list(np.linalg.norm(
        nodes[bars[i][1]] - nodes[bars[i][0]]) for i in range(len(bars)))   
    gL = 2 * len(nodes)
    mrg = np.zeros((gL, gL))  # matriz de rigidez global
    # calcula as matrizes de rigidez locais
    for ind, bar in enumerate(bars):
        # calcula o vetor de cossenos diretores
        C = (nodes[bar[1]] - nodes[bar[0]]) / comp_barras[ind]
        # matriz provisória para facilitar o indexing
        mrlb = np.zeros((4, 4))
        #mrl = np.zeros((4,4))
        # preenche o primeiro quadrante
        mrlb[0, 0] = C[0] * C[0]  # preenche a diagonal
        mrlb[1, 1] = C[1] * C[1]  # preenche a diagonal
        mrlb[1, 0] = C[0] * C[1]
        mrlb[0, 1] = C[0] * C[1]
        # preenche o segundo quadrante
        mrlb[2:4, 0:2] = -1 * mrlb[0:2, 0:2]
        # preenche o 3 e 4 quadrante
        mrlb[:, 2:4] = -1 * mrlb[:, 0:2]
        # print(mrlb)
        # multiplica pela area da barra
        # print("This is the bar {} and \n{}".format(ind+1, mrlb))
        mrlb = mrlb * area[ind]
        mrlb = mrlb/comp_barras[ind]
        # matriz identidade
        # traz as forcas referentes aos nós da barra
        forca = forcas_index[bar[0]] + forcas_index[bar[1]]
        # keep the identity of the forces
        mat_id = np.zeros((gL, 4))
        mat_id[forca, np.arange(4)] = 1
        prev = np.matmul(mat_id, mrlb)
        mrg += np.matmul(prev, mat_id.T)
        mrl[:, :, ind] = np.copy(mrlb)
        # multiplica pelo momento de inercia do material
    mrg = mrg * mom_e
    mrl = mrl * mom_e
    return (mrl, mrg)


def desloc(desloc_rest, bars, nodes, area, mom_e, loads):
    # matriz dos deslocamentos
    desloc = np.zeros((2*len(nodes), 1))
    # matriz das forças
    forcas = np.zeros((2*len(nodes), 1))
    # fills the loads
    for x, y in enumerate(loads):
        forcas[y[0]] = y[1]
    desloc_unrest = list(i for i in range(len(forcas)) if i not in desloc_rest)
    # calculates the local matrices
    mrl, mrg = stiff_matrices(bars, nodes, area, mom_e)
    du_xx, du_yy = np.meshgrid(desloc_unrest, desloc_unrest)
    # matriz de rigidez global nas coordenadas com desloc irrestrito
    k_du = mrg[du_xx, du_yy]
    # resolve os deslocamentos irrestritos
    desloc[desloc_unrest] = np.linalg.solve(k_du, forcas[desloc_unrest])
    return desloc


def bar_force(desloc_rest, bars, nodes, area, mom_e, loads):
    # matriz das forças
    forcas = np.zeros((2*len(nodes), 1))
    # fills the loads
    for x, y in enumerate(loads):
        forcas[y[0]] = y[1]
     # relaciona os nos as forcas
    forcas_index = []
    fi = 0
    for ff in range(len(nodes)):
        forcas_index.append((fi, fi+1))
        fi += 2
    desl = desloc(desloc_rest, bars, nodes, area, mom_e, loads)
    # calculates the global and local stiffness matrix
    mrl, mrg = stiff_matrices(bars, nodes, area, mom_e)
    # calcula as forcas nas barras
    forca_barra = np.zeros((len(bars), 4))
    for bb in range(len(bars)):
        # array provisório para facilitar o indexing
        mrlf = mrl[:, :, bb]
        # traz as coordenadas das forcas
        forca = list(forcas_index[bars[bb][0]] + forcas_index[bars[bb][1]])
        # multiplica a matriz local pelos deslocamentos
        forca_barra[bb] = np.reshape(np.matmul(mrlf, desl[forca]), 4)
    return forca_barra


def support_reac(desloc_rest, bars, nodes, area, mom_e, loads):
    disp = desloc(desloc_rest, bars, nodes, area, mom_e, loads)
    mrl, mrg = stiff_matrices(bars, nodes, area, mom_e)
    desloc_unrest = list(i for i in range(
        2 * len(nodes)) if i not in desloc_rest)
    # chama o setor correto
    ra_xx, ra_yy = np.meshgrid(desloc_unrest, desloc_rest)
    # corta a matriz de rigidez global
    ra_du = mrg[ra_xx, ra_yy]
    # calcula areacao dos apoios
    reacao_apoios = np.matmul(ra_du, disp[desloc_unrest])
    return reacao_apoios


def scalar_bar_force(desloc_rest, bars, nodes, area, mom_e, loads):
    forca_barra = bar_force(desloc_rest, bars, nodes, area, mom_e, loads)
    # retorna as forcas escalares
    forca_escalar = np.apply_along_axis(np.linalg.norm, 1, forca_barra[:, :2])
    return forca_escalar

def main():
    # Momento de inércia
    mom = 68950000.0
    # Nós da treliça
    nos = np.array([[0.0,0.0], [0.0, 9.144], [9.144, 0.0],[9.144,9.144],[18.288,0.0],[18.288,9.144]])
    # Barras da treliça
    barras = [[0,2],[0,3],[1,2],[1,3],[3,2],[3,4],[5,2],[3,5],[2,4],[4,5]]
    # Areas das secoes
    areas = [19519.11, 64.5, 14784.81, 9705.02, 65.48,
                               300.48, 4819.16, 13610.37, 13930.3, 64.5]
    a = stiff_matrices(barras, nos, areas, mom)
    # print(a[1])

if __name__ == '__main__':
    main()