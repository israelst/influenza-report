#coding:utf8

__author__ = 'Marcelo Ferreira da Costa Gomes'
from numpy import *
import scipy as sp
from pandas import *
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
from argparse import RawDescriptionHelpFormatter
import matplotlib.font_manager as fm
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import matplotlib.ticker as ticker


# Load R MEM package:
mem = importr('mem')
ro.r.require('mem')
# UF codes
tabela_ufnome = {11: 'Rondônia',
                 12: 'Acre',
                 13: 'Amazonas',
                 14: 'Roraima',
                 15: 'Pará',
                 16: 'Amapá',
                 17: 'Tocantins',
                 21: 'Maranhão',
                 22: 'Piauí',
                 23: 'Ceará',
                 24: 'Rio Grande do Norte',
                 25: 'Paraíba',
                 26: 'Pernambuco',
                 27: 'Alagoas',
                 28: 'Sergipe',
                 29: 'Bahia',
                 31: 'Minas Gerais',
                 32: 'Espírito Santo',
                 33: 'Rio de Janeiro',
                 35: 'São Paulo',
                 41: 'Paraná',
                 42: 'Santa Catarina',
                 43: 'Rio Grande do Sul',
                 50: 'Mato Grosso do Sul',
                 51: 'Mato Grosso',
                 52: 'Goiás',
                 53: 'Distrito Federal',
                 'AfAmN': 'Região 1',
                 'Aw': 'Região 3',
                 'AsBSh': 'Região 2',
                 'Cf': 'Região 4'}
tabela_ufcod = {v: k for k,v in tabela_ufnome.items()}
fontproplgd = fm.FontProperties('Oswald')
fontproplgd.set_size(28)
fontproplbl = fm.FontProperties('Oswald')
fontproplbl.set_size(42)
fontproplblinset = fm.FontProperties('Oswald')
fontproplblinset.set_size(30)
fontpropticks = fontproplblinset.copy()
fontpropticks.set_size(24)
fontpropticksinset = fontpropticks.copy()
fontpropticksinset.set_size(20)

def applymem(df):
    rdf = pandas2ri.py2ri(df)
    seasons = sort(list(df.columns.drop(['UF','isoweek'])))[:-1]
    # Discard 2009 season if present:
    seasons = sorted(set(seasons).difference(['SRAG2009']))
    rseasons = ro.StrVector(seasons)
    ro.globalenv['df'] = rdf
    ro.globalenv['seasons'] = rseasons
    # # Method for obtaining typical time series evolution (default 2)
    # ro.globalenv['par.type.curve'] = 2
    # # Method for obtaining pre/post-epidemic threshold (default 4)
    # ro.globalenv['par.type.threshold'] = 2
    # # Method for obtaining intensity thresholds (default 4)
    # ro.globalenv['par.type.intensity'] = 2
    # # Method for obtaining outbreak start and length (default 6)
    # ro.globalenv['par.type.other'] = 2
    # # Total number of points to obtain pre/post-threshold (will take n/seasons from each)
    # ro.globalenv['par.n.max'] = 30
    # # Confidence interval for modelled curve
    # ro.globalenv['par.level.curve'] = 0.90
    # # Confidence interval for pre/post-thresold
    # ro.globalenv['par.level.threshold'] = 0.95
    # # Quantiles for intensity thresholds
    # ro.globalenv['par.level.intensity'] = ro.FloatVector([0.40, 0.90, 0.975])
    #
    # epimemrslt = ro.r('memmodel(i.data=subset(df, select=seasons), i.type.curve=par.type.curve,' +
    #                   'i.type.threshold=par.type.threshold, i.type.intensity=par.type.intensity,' +
    #                   'i.type.other=par.type.other, i.n.max=par.n.max, i.level.curve=par.level.curve,' +
    #                   'i.level.threshold=par.level.threshold, i.level.intensity=par.level.intensity)')

    ro.globalenv['df'] = rdf
    ro.globalenv['seasons'] = rseasons
    ro.globalenv['par.type.curve'] = 2
    ro.globalenv['par.n.max'] = 20
    ro.globalenv['par.level.curve'] = 0.90
    ro.globalenv['par.level.threshold'] = 0.90

    epimemrslt = ro.r('memmodel(i.data=subset(df, select=seasons), i.type.curve=par.type.curve,' +
                      'i.n.max=par.n.max, i.level.curve=par.level.curve, i.level.threshold=par.level.threshold)')

    # Pre-epidemic threshold:
    epithreshold = pandas2ri.ri2py_dataframe(epimemrslt.rx2('pre.post.intervals')).loc[0,2]
    typrealcurve = pandas2ri.ri2py_dataframe(epimemrslt.rx2('typ.real.curve'))

    # Check for seasons below threshold:
    dropseasons = set()
    for s in seasons:
        if df[s].max() < epithreshold:
            dropseasons.add(s)
    # Drop seasons below threshold and rerun algorithm:
    episeasons = list(seasons)
    if len(dropseasons) > 0 and len(dropseasons) < len(seasons):
        episeasons = sorted(list(set(seasons).difference(dropseasons)))
        ro.globalenv['episeasons'] = ro.StrVector(episeasons)

        # epimemrslt = ro.r('memmodel(i.data=subset(df, select=episeasons), i.type.curve=par.type.curve,' +
        #                   'i.type.threshold=par.type.threshold, i.type.intensity=par.type.intensity,' +
        #                   'i.type.other=par.type.other, i.n.max=par.n.max, i.level.curve=par.level.curve,' +
        #                   'i.level.threshold=par.level.threshold, i.level.intensity=par.level.intensity)')

        epimemrslt = ro.r('memmodel(i.data=subset(df, select=episeasons), i.type.curve=par.type.curve,' +
                          'i.n.max=par.n.max, i.level.curve=par.level.curve, i.level.threshold=par.level.threshold)')

    # Store results in python dictionary of objects
    pyepimemrslt = {}
    rovector = [ro.vectors.StrVector, ro.vectors.IntVector, ro.vectors.FloatVector, ro.vectors.Vector]
    for name in epimemrslt.names:
        rdata = epimemrslt.rx2(name)
        if name == 'call':
            pyepimemrslt.update({name: str(rdata)})
        elif type(rdata) in rovector:
            pyepimemrslt.update({name: pandas2ri.ri2py_vector(rdata)})
        else:
            pyepimemrslt.update({name: pandas2ri.ri2py_dataframe(rdata)})

    # typ.curve is the typical curve obtained from averaging over epimemic seasons with time rescaled
    # so that the start of the epidemic period coincides with mean.start
    pyepimemrslt['typ.curve'].rename(columns={0:'baixo', 1:'mediano', 2:'alto'}, inplace=True)
    pyepimemrslt['typ.curve']['mediano'].fillna(0, inplace=True)
    pyepimemrslt['typ.curve']['baixo'] = pyepimemrslt['typ.curve']['baixo'].where(pyepimemrslt['typ.curve']['baixo']>=0,
                                                                                  other=0)
    pyepimemrslt['typ.curve']['baixo'] = pyepimemrslt['typ.curve']['baixo'].\
        where( (-pyepimemrslt['typ.curve']['baixo'].isnull()), other=pyepimemrslt['typ.curve']['mediano'])
    pyepimemrslt['typ.curve']['alto'] = pyepimemrslt['typ.curve']['alto'].\
        where((-pyepimemrslt['typ.curve']['alto'].isnull()), other=pyepimemrslt['typ.curve']['mediano'])
    pyepimemrslt['pre.post.intervals'].rename(index={0: 'pre', 1: 'post'}, inplace=True)

    # typ.real.curve is the typical curve without time shift, that is, respecting the original weeks from data
    # this curve is better to keep all seasons, not only the epidemic ones.
    pyepimemrslt['typ.real.curve'] = typrealcurve.copy()
    pyepimemrslt['typ.real.curve'].rename(columns={0:'baixo', 1:'mediano', 2:'alto'}, inplace=True)
    pyepimemrslt['typ.real.curve']['mediano'].fillna(0, inplace=True)
    pyepimemrslt['typ.real.curve']['baixo'] = pyepimemrslt['typ.real.curve']['baixo'].\
        where(pyepimemrslt['typ.real.curve']['baixo']>=0, other=0)
    pyepimemrslt['typ.real.curve']['baixo'] = pyepimemrslt['typ.real.curve']['baixo'].\
        where( (-pyepimemrslt['typ.real.curve']['baixo'].isnull()), other=pyepimemrslt['typ.real.curve']['mediano'])
    pyepimemrslt['typ.real.curve']['alto'] = pyepimemrslt['typ.real.curve']['alto'].\
        where((-pyepimemrslt['typ.real.curve']['alto'].isnull()), other=pyepimemrslt['typ.real.curve']['mediano'])
    newcols = {}
    for k,v in enumerate(episeasons):
        newcols[k] = str(v) + ' transladado'
    pyepimemrslt['moving.epidemics'].rename(columns=newcols, inplace=True)

    return pyepimemrslt

def main(fname, sep=',', uflist='all'):

    df = pd.read_csv(fname, sep=sep)
    dfinset = pd.read_csv(fname.replace('-incidence',''), sep=sep)
    if 'Região' in list(df.columns):
        df.rename(columns={'Região': 'UF'}, inplace=True)
        dfinset.rename(columns={'Região': 'UF'}, inplace=True)

    plt.interactive(False)
    if uflist == 'all':
        uflist = list(df.UF.unique())

    for uf in uflist:
        try:
            uf = int(uf)
        except:
            pass
        if uf not in list(df.UF.unique()):
            continue
        print(uf)
        dftmp = df[df.UF == uf].reset_index().drop('index', axis=1).copy()
        dftmpinset = dfinset[dfinset.UF == uf].reset_index().drop('index', axis=1).copy()
        seasons = sort(list(dftmp.columns.drop(['UF','isoweek'])))
        lastseason = seasons[-1]
        seasons = np.delete(seasons, -1)
        sns.set_style('darkgrid')
        sns.set_context("talk")
        sns.set_palette('Set2', len(seasons)+4)
        colorcode = sns.color_palette('Set2', len(seasons)+4)

        try:
            thresholds = applymem(dftmp)
            thresholdsinset = applymem(dftmpinset)
            dftmp['mediana pré-epidêmica'] = thresholds['pre.post.intervals'].loc['pre',1]
            dftmp['limiar pré-epidêmico'] = thresholds['pre.post.intervals'].loc['pre',2]
            dftmp['limiar pós-epidêmico'] = thresholds['pre.post.intervals'].loc['post',2]
            dftmp['intensidade baixa'] = thresholds['epi.intervals'].loc[0,3]
            dftmp['intensidade alta'] = thresholds['epi.intervals'].loc[1,3]
            dftmp['intensidade muito alta'] = thresholds['epi.intervals'].loc[2,3]
            dftmp['corredor baixo'] = thresholds['typ.real.curve']['baixo']
            dftmp['corredor mediano'] = thresholds['typ.real.curve']['mediano']
            dftmp['corredor alto'] = thresholds['typ.real.curve']['alto']
            dftmp['se relativa ao início do surto'] = dftmp['isoweek'] - thresholds['mean.start'][0]
            dftmp['se típica do início do surto'] = thresholds['mean.start'][0]
            dftmp['duração típica do surto'] = thresholds['mean.length'][0]
            dftmp['curva epi. baixa'] = thresholds['typ.curve']['baixo']
            dftmp['curva epi. mediana'] = thresholds['typ.curve']['mediano']
            dftmp['curva epi. alta'] = thresholds['typ.curve']['alto']
            epicols = list(thresholds['moving.epidemics'].columns)
            dftmp[epicols] = thresholds['moving.epidemics']
            dftmp['n.seasons'] = thresholds['n.seasons'][0]
            dftmp.to_csv('%s-mem-incidencia.csv' % tabela_ufnome[uf].replace(' ','_'), index=False)

            dftmpinset['limiar pré-epidêmico'] = thresholdsinset['pre.post.intervals'].loc['pre',2]
            dftmpinset['limiar pós-epidêmico'] = thresholdsinset['pre.post.intervals'].loc['post',2]
            dftmpinset['intensidade baixa'] = thresholdsinset['epi.intervals'].loc[0,3]
            dftmpinset['intensidade alta'] = thresholdsinset['epi.intervals'].loc[1,3]
            dftmpinset['intensidade muito alta'] = thresholdsinset['epi.intervals'].loc[2,3]
            dftmpinset['corredor baixo'] = thresholdsinset['typ.real.curve']['baixo']
            dftmpinset['corredor mediano'] = thresholdsinset['typ.real.curve']['mediano']
            dftmpinset['corredor alto'] = thresholdsinset['typ.real.curve']['alto']
            dftmpinset['se relativa ao início do surto'] = dftmpinset['isoweek'] - thresholdsinset['mean.start'][0]
            dftmpinset['se típica do início do surto'] = thresholdsinset['mean.start'][0]
            dftmpinset['duração típica do surto'] = thresholdsinset['mean.length'][0]
            dftmpinset['curva epi. baixa'] = thresholdsinset['typ.curve']['baixo']
            dftmpinset['curva epi. mediana'] = thresholdsinset['typ.curve']['mediano']
            dftmpinset['curva epi. alta'] = thresholdsinset['typ.curve']['alto']
            epicols = list(thresholdsinset['moving.epidemics'].columns)
            dftmpinset[epicols] = thresholdsinset['moving.epidemics']
            dftmpinset['n.seasons'] = thresholdsinset['n.seasons'][0]
            dftmpinset.to_csv('%s-mem.csv' % tabela_ufnome[uf].replace(' ','_'), index=False)
            print(uf)


            fig, ax = plt.subplots(nrows=2, ncols=1, figsize = [20, 20])
            plt.subplots_adjust(hspace=0.3)

            # Set ymax at least = 1:
            if dftmp[list(set(seasons).union([lastseason]))].max().max() < 1 or uf in [32,33]:
                ax[0].set_ylim([0,1])
                ax[1].set_ylim([0,1])
            # if uf == 33:
            #     ax[0].set_ylim([0,0.25])
            # elif uf == 32:
            #     ax[0].set_ylim([0,0.3])

            ax[0].fill_between(dftmp['isoweek'], 0, dftmp['corredor baixo'], color='green', alpha=0.5)
            ax[0].fill_between(dftmp['isoweek'], dftmp['corredor baixo'], dftmp['corredor mediano'], color='yellow',
                               alpha=0.5)
            ax[0].fill_between(dftmp['isoweek'], dftmp['corredor mediano'], dftmp['corredor alto'], color='orange',
                               alpha=0.5)

            dftmp.plot(ax=ax[0], x='isoweek', y=seasons)
            dftmp.plot(ax=ax[0], x='isoweek', y=lastseason, color='k', lw=3)
            dftmp.plot(ax=ax[0], x='isoweek', y='limiar pré-epidêmico', style='--', color='red', alpha=0.8)
            dftmp.plot(ax=ax[0], x='isoweek', y='intensidade baixa', style='--')
            dftmp.plot(ax=ax[0], x='isoweek', y='intensidade alta', style='--')
            dftmp.plot(ax=ax[0], x='isoweek', y='intensidade muito alta', style='--', color=colorcode[-1])

            # Check for maximum value on y-axis and fill from 'corredor alto' to maxy
            dftmp.plot(ax=ax[0], x='isoweek', y='corredor alto', legend=False, alpha=0)
            miny, maxy = ax[0].get_ylim()
            del(ax[0].lines[-1])
            ax[0].fill_between(dftmp['isoweek'], dftmp['corredor alto'], maxy, color='red', alpha=0.5)
            ax[0].set_ylim([miny, maxy])
            for label in ax[0].get_xticklabels() :
                label.set_fontproperties(fontpropticks)
            for label in ax[0].get_yticklabels() :
                label.set_fontproperties(fontpropticks)

            #### Start absolute value plot as inset ####
            axinset = inset_axes(ax[0], width='35%', height='35%', loc=1)
            if uf == 33:
                axinset.set_ylim([0,40])
            elif uf == 32:
                axinset.set_ylim([0,12])
            dftmpinset.plot(ax=axinset, x='isoweek', y=seasons)
            dftmpinset.plot(ax=axinset, x='isoweek', y=lastseason, color='k', lw=3)
            dftmpinset.plot(ax=axinset, x='isoweek', y='limiar pré-epidêmico', style='--', color='red', alpha=0.8)
            axinset.legend_.remove()
            axinset.set_xlabel('SE', fontproperties = fontproplblinset)
            axinset.set_ylabel('Casos', fontproperties = fontproplblinset)
            axinset.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
            for label in axinset.get_xticklabels() :
                label.set_fontproperties(fontpropticksinset)
            for label in axinset.get_yticklabels() :
                label.set_fontproperties(fontpropticksinset)


            #### Start plot relative to outbreak typical curve ####

            ax[1].fill_between(dftmp['se relativa ao início do surto'], 0, dftmp['curva epi. baixa'], color='green',
                               alpha=0.5)
            ax[1].fill_between(dftmp['se relativa ao início do surto'], dftmp['curva epi. baixa'],
                               dftmp['curva epi. mediana'], color='yellow', alpha=0.5)
            ax[1].fill_between(dftmp['se relativa ao início do surto'], dftmp['curva epi. mediana'],
                               dftmp['curva epi. alta'], color='orange', alpha=0.5)
            dftmp.plot(ax=ax[1], x='se relativa ao início do surto', y='curva epi. mediana', color='silver',
                       label='tendência mediana')
            dftmp.plot(ax=ax[1], x='se relativa ao início do surto', y='limiar pré-epidêmico', style='--',
                       color='red', alpha=0.8)
            dftmp.plot(ax=ax[1], x='se relativa ao início do surto', y='limiar pós-epidêmico', style='--',
                       color='green', alpha=0.5)

            epicolor = []
            for s in epicols:
                s = s.strip(' transladado')
                n = list(seasons).index(s)
                epicolor.append(colorcode[n])
            dftmp.plot(ax=ax[1], x='se relativa ao início do surto', y=epicols, color=epicolor)
            # Check for maximum value on y-axis and fill from 'corredor alto' to maxy
            dftmp.plot(ax=ax[1], x='se relativa ao início do surto', y='curva epi. alta', legend=False, alpha=0)
            miny, maxy = ax[1].get_ylim()
            del(ax[1].lines[-1])
            ax[1].fill_between(dftmp['se relativa ao início do surto'], dftmp['curva epi. alta'], maxy, color='red',
                               alpha=0.5)
            ax[1].set_ylim([miny, maxy])
            ax[1].plot([0,0], [miny, maxy], '--', color='silver')
            duracao = int(thresholds['mean.length'][0])
            ax[1].plot([duracao, duracao], [miny, maxy], '--', color='silver')

            ax[1].set_title('Tendência ao longo do surto', fontproperties=fontproplbl)
            epistart = int(thresholds['mean.start'][0])
            ax[1].set_xlabel('SE em relação à semana típica de início do surto (SE=%s)' % epistart,
                             fontproperties=fontproplbl)
            minx, maxx = ax[1].get_xlim()
            xticks = sort(np.append(np.arange(0,int(minx),-4),np.arange(4,int(maxx),4)))
            ax[1].set_xticks(xticks)
            ax[1].set_xticklabels(xticks, fontproperties=fontpropticks)
            for label in ax[0].get_yticklabels() :
                label.set_fontproperties(fontpropticks)
            ax[1].set_ylabel('Incidência (por 100mil habitantes)', fontproperties=fontproplbl)
            box = ax[1].get_position()
            ax[1].set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax[1].legend(prop=fontproplgd, loc='center left', bbox_to_anchor=(1,0.5))

        except:
            dftmp.to_csv('%s-memfailed-incidencia.csv' % tabela_ufnome[uf].replace(' ','_'), index=False)
            fig, axi = plt.subplots(nrows=1, ncols=1, figsize = [20, 10])
            ax = [axi]

            if dftmp[list(set(seasons).union([lastseason]))].max().max() < 1 or uf in [32,33]:
                ax[0].set_ylim([0,1])

            dftmp.plot(ax=ax[0], x='isoweek', y=seasons)
            dftmp.plot(ax=ax[0], x='isoweek', y=lastseason, color='k', lw=3)

            #### Start absolute value plot as inset ####
            sns.set_style('whitegrid')
            axinset = inset_axes(ax[0], width='35%', height='35%', loc=1)
            dftmpinset.plot(ax=axinset, x='isoweek', y=seasons)
            dftmpinset.plot(ax=axinset, x='isoweek', y=lastseason, color='k', lw=3)
            axinset.legend_.remove()
            axinset.set_xlabel('SE')
            axinset.set_ylabel('Casos')
            axinset.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

        ax[0].set_title(tabela_ufnome[uf], fontproperties=fontproplbl)
        ax[0].set_xlabel('SE', fontproperties=fontproplbl)
        ax[0].set_ylabel('Incidência (por 100mil habitantes)', fontproperties=fontproplbl)
        xticks = np.arange(4,53,4)
        ax[0].set_xticks(xticks)
        ax[0].set_xticklabels(xticks)
        # Shrink current axis by 10%
        box = ax[0].get_position()
        ax[0].set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax[0].legend(prop=fontproplgd, loc='center left', bbox_to_anchor=(1,0.5))

        #plt.tight_layout()
        plt.savefig(tabela_ufnome[uf].replace(' ','_')+'-inset.svg')
        plt.clf()
        plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate MEM analysis from cleaned SINAN-SRAG data,\n" +
                                     "for specified Federal Units, if any. If none specified, runs for all.\n" +
                                     "Example usage:\n" +
                                     "python3 sinan-mem-inset.py --path clean_data4mem-regiao-incidence.csv --uflist" +
                                     " Aw Cf\n",
                                     formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--path', help='Path to data file')
    parser.add_argument('--sep', help='Column separator', default=',')
    parser.add_argument('--uflist', nargs='*', default='all')
    args = parser.parse_args()
    main(args.path, args.sep, args.uflist)
