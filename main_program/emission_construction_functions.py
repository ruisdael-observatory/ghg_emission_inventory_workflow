
"""
Functions used to create point source and area emissions:

"""


def determine_snap(subdoelgroep):
    
    if subdoelgroep in ['Opwekking electriciteit']:
        
        isnap = 1
        
    elif subdoelgroep in ['Energiegebruik en processen Handel, Diensten en Overheid (HDO)',
                          'Energiegebruik en processen Handel. Diensten en Overheid (HDO)',
                          'Indirecte emissies broeikasgassen Consumenten',
                          'Productgebruik Consumenten',
                          'Energiegebruik Consumenten',
                          'Indirecte emissies broeikasgassen Handel, Diensten en Overheid (HDO)']:
        isnap = 2
        
    elif subdoelgroep in ['Overig bouw',
                          'Chemische Industrie basisproducten',
                          'Metaalelektro',
                          'Industrie overig',
                          'Indirecte emissies broeikasgassen Bouw',
                          'Rubber- en kunstofverw. Industrie',
                          'Chemische Industrie overig',
                          'Chemische Industrie kunstmeststoffen',
                          'Houtbewerkende industrie',
                          'Bouwmaterialenindustrie',
                          'Voedings- en genotmiddelenindustrie',
                          'Lederindustrie',
                          'Grafische industrie',
                          'Papier(waren)',
                          'Textiel- en tapijtindustrie',
                          'Basismetaal',
                          'Drinkwaterbedrijven',
                          'Chemische Industrie bestrijdingsmiddelen',
                          'Overige processen']:
        isnap = 3
        
    elif subdoelgroep in ['Raffinage en verwerking',
                          'Railverkeer', #note:Railverkeer is snap 8, but since we refine it with shipping lines Railverkeer is put here as snap 4 has no refinement
                          'Mobiele werktuigen']: #note:Mobiele werktuigen is snap 8, but since we refine it with shipping lines Mobiele werktuigen is put here as snap 4 has no refinement
        
        isnap = 4
    
    elif subdoelgroep in ['Transport en distributie olie en gas',
                          'Olie- gaswinning land',
                          'Olie- gaswinning continentaal plat']:
        isnap = 5
        
    elif subdoelgroep in ['Smeermiddelengebruik-verkeer',
                          'Wegverkeer - uitlaatgassen',
                          'Wegverkeer - niet uitlaatgassen']:
        isnap = 7
        
    elif subdoelgroep in [
                          'Binnenscheepvaart',
                          'Zeescheepvaart NCP (inclusief ankerliggers)',
                          'Zeescheepvaart stilliggend',
                          'Zeescheepvaart varend op Nederlands grondgebied',
                          'Recreatievaart',
                          'Visserij']:
        isnap = 8 #only shipping for now
        
    elif subdoelgroep in ['Indirecte emissies broeikasgassen Afvalverwijdering',
                          "AVI's",
                          'Energiegebruik en processen Riolering en waterzuiveringsinstallaties',
                          'Overige afvalbedrijven',
                          'Storten',
                          'Huishoudelijk afvalwater',
                          'Indirecte emissies broeikasgassen Riolering en waterzuiveringsinstallaties']:
        isnap = 9
        
    elif subdoelgroep in ['Indirecte emissies broeikasgassen Landbouw',
                          'Productgebruik Landbouw',
                          'Energiegebruik Landbouw',
                          'Landbouwbedrijven',
                          'Kunstmest',
                          'Landbouwhuisdieren - pluimvee',
                          'Landbouwhuisdieren - overige dieren',
                          'Landbouwhuisdieren - rundvee',
                          'Bodems natuur',
                          'Bodems landbouw',
                          'Landbouwhuisdieren - varkens',
                          'Overig drinkwater',
                          'Composteren',
                          'Particuliere landbouwactiviteiten']:
        isnap = 10
    
    elif subdoelgroep in ['Luchtvaart']:
        isnap = 11
    else:
        print(f'Unknown subdoelgroep: {subdoelgroep}')
        isnap = 0
        
    return isnap

def emissieoorzaak_snap(emissieoorzaak):
       
    emissieoorzaak_snap1 = ['SBI 35.111 (per bedrijf): Elektriciteitsproduktie',
                            'SBI 35 (per bedrijf): Productie en distributie van elektriciteit en gas']
    
    emissieoorzaak_snap2 = [
 'SBI 46/47 (per bedrijf): Detail- en groothandel',
 'SBI 49-53 (per bedrijf): Transport, communicatie',
 'SBI 64 (per bedrijf): Financiele dienstverlening (excl. verzekeringen en pensioenfondsen)',
 'SBI 64.9 (per bedrijf): Overige financiele dienstverlening',
 'SBI 85-88 (per bedrijf): Onderwijs en gezondheids- en welzijnszorg', 
 'SBI 59/60/90/91/93/96 (per bedrijf): Cultuur, sport, recreatie en overige dienstverlening', 
 'SBI 18/58 (per bedrijf): Uitgeverijen, drukkerijen, reproductie van opgenomen media', 
 'SBI 63 (per bedrijf) - Gegevensverwerking, webhosting en aanverwante activiteiten', 
 'SBI 84.1 (per bedrijf): Openbaar bestuur']
    
    
    emissieoorzaak_snap3 = [
 'SBI 10-12 (per bedrijf): Voedings- & genotmiddelenindustrie',
 'SBI 10.1 (per bedrijf): Slachterijen en vleesverwerking',
 'SBI 10.3 (per bedrijf): Groente- en fruitverwerking',
 'SBI 10.5 (per bedrijf): Zuivelindustrie',
 'SBI 10.6 (per bedrijf): Meelproduktie (excl. zetmeel)',
 'SBI 10.8 (per bedrijf): Overige voedingsmiddelenindustrie (exclusief SBI 10.81 en 10.82)',
 'SBI 10.9 (per bedrijf): Diervoederindustrie',
 'SBI 11.07 (per bedrijf): Vervaardiging van dranken',
 'SBI 13/14 (per bedrijf): Vervaardiging van textiel en kleding',
 'SBI 15.11 (per bedrijf): Leerlooierijen en bontbereiding',
 'SBI 16 (per bedrijf): Houtindustrie en vervaardiging van artikelen van hout, kurk, riet en vlechtwerk (geen meubels)',
 'SBI 17 (per bedrijf): Vervaardiging van papier, karton en papier- en kartonwaren',
 'SBI 17.1 (per bedrijf): Vervaardiging van papierpulp, papier en karton',
 'SBI 17.2 (per bedrijf): Vervaardiging van papier- en kartonwaren',
 'SBI 20.1 (per bedrijf): Vervaardiging van chemische basisproducten',
 'SBI 20.11 (per bedrijf): Vervaardiging van industriÃ«le gassen',
 'SBI 20.11 (per bedrijf): Vervaardiging van industriële gassen',
 'SBI 20.11 (per bedrijf): Vervaardiging van industriÃle gassen',
 'SBI 20.11 (per bedrijf): Vervaardiging van industriÎle gassen',      
 'SBI 20.13 (per bedrijf): Basischemie anorganisch',
 'SBI 20.141 (per bedrijf): Vervaardiging van petrochemische producten',
 'SBI 20.149 (per bedrijf): Basischemie organisch (geen petrochemische producten)',
 'SBI 20.15 (per bedrijf): Vervaardiging van kunstmeststoffen en stikstofverbindingen',
 'SBI 20.16 (per bedrijf): Vervaardiging van kunststof in primaire vorm',
 'SBI 20.2 (per bedrijf): Chemische bestrijdingsmiddelenindustrie',
 'SBI 20.3 (per bedrijf): Vervaardiging van verf, lak, vernis, inkt en mastiek',
 'SBI 20.4 (per bedrijf): Vervaardiging was- en schoonmaakmiddelen, parfums en cosmetica',
 'SBI 20.5 (per bedrijf): Overige chemische producten',
 'SBI 20.59 (per bedrijf): Vervaardiging van overige chemische producten n.e.g.',
 'SBI 20.6 (per bedrijf): Vervaardiging van synthetische en kunstmatige vezels',
 'SBI 21.1 (per bedrijf): Vervaardiging van farmaceutische producten',
 'SBI 21.20 (per bedrijf): Vervaardiging van farmaceutische producten (geen grondstoffen)',
 'SBI 22.1 (per bedrijf): Vervaardiging van producten van rubber',
 'SBI 22.2 (per bedrijf): Vervaardiging van producten van kunststof',
 'SBI 23 (per bedrijf): Bouwmaterialen- en glasindustrie',
 'SBI 23.1 (per bedrijf): Vervaardiging van glas en glaswerk',
 'SBI 23.2-23.4 (per bedrijf): Vervaardiging van keramische producten',
 'SBI 23.32 (per bedrijf): Vervaardiging van bakstenen en dakpannen',
 'SBI 23.51 (per bedrijf): Vervaardiging van cement',
 'SBI 23.9 (per bedrijf): Vervaardiging van overige niet-metaalhoudende minerale producten',
 'SBI 24 (per bedrijf): Vervaardiging van metalen in primaire vorm',
 'SBI 24.1-24.3 (per bedrijf): Basismetaalindustrie, verwerking en vervaardiging ijzer en staal',
 'SBI 24.4/24.5 (per bedrijf): Basismetaalindustrie, vervaardiging van non-ferro metalen en gieten van metalen',
 'SBI 24.45 (per bedrijf) Vervaardiging van overige non-ferrometalen, aluminium',
 'SBI 24.45 (per bedrijf) Vervaardiging van overige non-ferrometalen, koper',
 'SBI 24.45 (per bedrijf) Vervaardiging van overige non-ferrometalen, lood',
 'SBI 24.45 (per bedrijf) Vervaardiging van overige non-ferrometalen, zink',
 'SBI 24.5 (per bedrijf): Gieten van metalen',
 'SBI 25 (per bedrijf): Metaalproductenindustrie (exclusief machinebouw)',
 'SBI 26/27 (per bedrijf): Elektrotechnische industrie',
 'SBI 28 (per bedrijf): Machinebouw',
 'SBI 29 (per bedrijf): Auto-industrie',
 'SBI 30 (per bedrijf): Overige transportmiddelen',
 'SBI 30.1 (per bedrijf): Scheepsbouw',
 'SBI 31/32 (per bedrijf): Vervaardiging van meubels en overige goederen',
 "SBI 45 (per bedrijf): Handel en reparatie van auto's en motorfietsen",  
 'SBI 20.12 (per bedrijf): Vervaardiging van kleur- en verfstoffen', 
 'SBI 24.2 (per bedrijf): Vervaardiging van stalen buizen en pijpen',
 'SBI 20.41 (per bedrijf): Vervaardiging van was- en schoonmaakmiddelen',
 'SBI 24.3 (per bedrijf): Overige eerste verwerking van ijzer en staal']
    
    emissieoorzaak_snap4 = ['SBI 19.201 (per bedrijf): Aardolieraffinage', 'SBI 15 (per bedrijf): Lederindustrie en bontbereiding', 
                            #'SBI 10.4 (per bedrijf): produktie oliÃn en vetten'
                           'SBI 10.4 (per bedrijf): produktie oliÃ«n en vetten',
                           'SBI 10.4 (per bedrijf): produktie oliÎn en vetten',
                           'SBI 10.4 (per bedrijf): produktie oliën en vetten',
                           'SBI 20.52 (per bedrijf): Vervaardiging van lijm en bereide kleefmiddelen']
    
    emissieoorzaak_snap5 = ['SBI 35.12 (per bedrijf): Transportnet voor aardgas', 
                            'SBI 08 (per bedrijf): Winning van delfstoffen (geen olie en gas)', 
                            'SBI 06/09.1 (per bedrijf): Aardolie- en gaswinning en dienstverlening voor de aardolie- en aardgaswinning']
    
    emissieoorzaak_snap8 = ['SBI 41-43 (per bedrijf): Bouwnijverheid',
                            'SBI 52.10/52.24 (per bedrijf): Laad-, los- en overslagactiviteiten en opslag']
    
    emissieoorzaak_snap9 = ['SBI 37: Afvalwaterinzameling en -behandeling',
                            "SBI 38.2 (per bedrijf): Afvalinzameling/beh, AVI's",
                            'SBI 38.2 (per bedrijf): Behandeling van afval',
                            'SBI 38.3 (per bedrijf): Voorbereiding tot recycling',
                            'RWZI spui van slibgistingsgas, CO2, individueel', 
                            'SBI 37 (per bedrijf): Afvalwaterinzameling en -behandeling']
    
    emissieoorzaak_snap10 = ['SBI 01 (per bedrijf): Landbouw, jacht en dienstverlening voor de landbouw en jacht']
    
    emissieoorzaak_snap11 = ['Vliegverkeer, Take Off', 
                            'Vliegverkeer, Climb Out',
                            'Vliegverkeer, Approach', 
                            'Vliegverkeer, Idle',
                            'Vliegverkeer, APU', 
                            'Vliegverkeer, GSE']
    


    if emissieoorzaak in emissieoorzaak_snap1:
        snap = 1
    elif emissieoorzaak in emissieoorzaak_snap2:
        snap = 2
    elif emissieoorzaak in emissieoorzaak_snap3:
        snap = 3
    elif emissieoorzaak in emissieoorzaak_snap4:
        snap = 4
    elif emissieoorzaak in emissieoorzaak_snap5:
        snap = 5
    elif emissieoorzaak in emissieoorzaak_snap8:
        snap = 8        
    elif emissieoorzaak in emissieoorzaak_snap9:
        snap = 9   
    elif emissieoorzaak in emissieoorzaak_snap10:
        snap = 10
    elif emissieoorzaak in emissieoorzaak_snap11:
        snap = 11
    else:
        snap = -1
        print('Unknown emissieoorzaak:', emissieoorzaak)
        
    return snap
    
    
    


def reademisoptions(rundir, show_log=False):
    
    from datetime import datetime
    from datetime import timedelta

    """
    Create variables required for the creation of hourly emission files, i.e.:
    1. Datetime object for iteration
    2. Domain specifications
    3. Tracer specifications

    Required input: textfile with combinations of keywords and values. Example:

    ----
    startyear   2017
    startmonth     1
    startday       1
    starthour      0
    runlength     24

    xmin 115000
    xmax 125000
    ymin 410000
    ymax 456700

    tracer co2 verkeer 1 power 7 residential 4
    tracer ch4 landbouw 10 waste 9
    ----
    """
    
    tracers, sources, categos = [], [], []

    with open(rundir, "r") as fh:
        for line in fh.readlines():
            words = line.split()

            if len(words) > 0:
                keyword = words[0]
                values  = words[1:]

                if keyword == 'startyear':
                    startyear  = int(values[0])
                elif keyword == 'startmonth':
                    startmonth = int(values[0])
                elif keyword == 'startday':
                    startday   = int(values[0])
                elif keyword == 'starthour':
                    starthour  = int(values[0])
                elif keyword == 'runlength':
                    runlength  = int(values[0])

                elif keyword == 'xmin':
                    xmin = float(values[0])
                elif keyword == 'xmax':
                    xmax = float(values[0])
                    
                elif keyword == 'ymin':
                    ymin = float(values[0])
                elif keyword == 'ymax':
                    ymax = float(values[0])
                    
                elif keyword == 'zmin':
                    zmin = int(values[0])
                elif keyword == 'zmax':
                    zmax = float(values[0])
                    
                elif keyword == 'tracer':
                    tracers.append(values[0])
                    sources.append(values[1::2])
                    categos.append(values[2::2])

    domainbounds = {'xmin': xmin, 'xmax': xmax, 
                    'ymin': ymin, 'ymax': ymax, 
                    'zmin': zmin, 'zmax': zmax}
    
    if show_log:
        for it, tracer in enumerate(tracers):
            print("Tracer: ", tracer, )
            for isrc, source in enumerate(sources[it]):
                print("{:02d}".format(int(categos[it][isrc])), source)
            print()
         
    tstart = datetime(startyear, startmonth, startday, starthour)
    tend   = tstart + timedelta(hours=runlength)

    if show_log:
        print("start:", tstart.strftime("%Y-%m-%d %H:%M"))
        print("end  :", tend  .strftime("%Y-%m-%d %H:%M"), "\n")

    return domainbounds, tstart, tend, tracers, sources, categos




def loadsnap():
    
    import numpy as np
    
    """
    Create arrays for temporal dissaggregation based on Denier van der Gon et al. (2010)
    using 10 SNAP cstegories
    
    INPUT: None
    OUTPUT: 
    
    tprof_mnth, np.array 10x12
    tprof_week, np.array 10x7
    tprof_hour, np.array 10x24
    
    """
    nsnap = 11

    tprof_mnth = np.zeros([nsnap, 12])
    tprof_week = np.zeros([nsnap, 7])
    tprof_hour = np.zeros([nsnap, 24])

    # SNAP 1 - Power generation
    tprof_mnth[0] = [1.20, 1.15, 1.05, 1.00, 0.90, 0.85, 0.80, 0.87, 0.95, 1.00, 1.08, 1.15]
    tprof_week[0] = [1.06, 1.06, 1.06, 1.06, 1.06, 0.85, 0.85]
    tprof_hour[0] = [0.79, 0.72, 0.72, 0.71, 0.74, 0.80, 0.92, 1.08, 1.19, 1.22, 1.21, 1.21, 
                     1.17, 1.15, 1.14, 1.13, 1.10, 1.07, 1.04, 1.02, 1.02, 1.01, 0.96, 0.88]

    # SNAP 2 - Residential, commercial and other combustion
    tprof_mnth[1] = [1.70, 1.50, 1.30, 1.00, 0.70, 0.40, 0.20, 0.40, 0.70, 1.05, 1.40, 1.65]
    tprof_week[1] = [1.08, 1.08, 1.08, 1.08, 1.08, 0.8, 0.8]
    tprof_hour[1] = [0.38, 0.36, 0.36, 0.36, 0.37, 0.50, 1.19, 1.53, 1.57, 1.56, 1.35, 1.16,
                     1.07, 1.06, 1.00, 0.98, 0.99, 1.12, 1.41, 1.52, 1.39, 1.35, 1.00, 0.42]

    # SNAP 3 - Industrial combustion
    tprof_mnth[2] = [1.10, 1.08, 1.05, 1.00, 0.95, 0.90, 0.93, 0.95, 0.97, 1.00, 1.02, 1.05]
    tprof_week[2] = [1.08, 1.08, 1.08, 1.08, 1.08, 0.8, 0.8]
    tprof_hour[2] = [0.75, 0.75, 0.78, 0.82, 0.88, 0.95, 1.02, 1.09, 1.16, 1.22, 1.28, 1.30,
                     1.22, 1.24, 1.25, 1.16, 1.08, 1.01, 0.95, 0.90, 0.85, 0.81, 0.78, 0.75]

    # SNAP 4 - Industrial processes
    tprof_mnth[3] = [1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.00, 0.84, 1.02, 1.02, 1.02, 0.90]
    tprof_week[3] = [1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.00]
    tprof_hour[3] = 24*[1.]

    # SNAP 5 - Extraction & distribution of fossil fuels
    tprof_mnth[4] = [1.20, 1.20, 1.20, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 1.20, 1.20, 1.20]
    tprof_week[4] =  7*[1.]
    tprof_hour[4] = 24*[1.]

    # SNAP 6 - Solvent use
    tprof_mnth[5] = [0.95, 0.96, 1.02, 1.00, 1.01, 1.03, 1.03, 1.01, 1.04, 1.03, 1.01, 0.91]
    tprof_week[5] = [1.2, 1.2, 1.2, 1.2, 1.2, 0.5, 0.5]
    tprof_hour[5] = [0.50, 0.35, 0.20, 0.10, 0.10, 0.20, 0.75, 1.25, 1.40, 1.50, 1.50, 1.50,
                     1.50, 1.50, 1.50, 1.50, 1.50, 1.40, 1.25, 1.10, 1.00, 0.90, 0.80, 0.70]

    # SNAP 7 - Road transport
    tprof_mnth[6] = [0.88, 0.92, 0.98, 1.03, 1.05, 1.06, 1.01, 1.02, 1.06, 1.05, 1.01, 0.93]
    tprof_week[6] = [1.02, 1.06, 1.08, 1.10, 1.14, 0.81, 0.79]
    tprof_hour[6] = [0.19, 0.09, 0.06, 0.05, 0.09, 0.22, 0.86, 1.84, 1.86, 1.41, 1.24, 1.20,
                     1.32, 1.44, 1.45, 1.59, 2.03, 2.08, 1.51, 1.06, 0.74, 0.62, 0.61, 0.44]

    # SNAP 8 - Other mobile sources and machinery
    tprof_mnth[7] = [0.88, 0.92, 0.98, 1.03, 1.05, 1.06, 1.01, 1.02, 1.06, 1.05, 1.01, 0.93]
    tprof_week[7] =  7*[1.]
    tprof_hour[7] = [0.19, 0.09, 0.06, 0.05, 0.09, 0.22, 0.86, 1.84, 1.86, 1.41, 1.24, 1.20,
                     1.32, 1.44, 1.45, 1.59, 2.03, 2.08, 1.51, 1.06, 0.74, 0.62, 0.61, 0.44]

    # SNAP 9 - Waste treatment and disposal
    tprof_mnth[8] = 12*[1.]
    tprof_week[8] =  7*[1.]
    tprof_hour[8] = 24*[1.]

    # SNAP 10 - Agriculture
    tprof_mnth[9] = [0.45, 1.30, 2.35, 1.70, 0.85, 0.85, 0.85, 1.00, 1.10, 0.65, 0.45, 0.45]
    tprof_week[9] =  7*[1.]
    tprof_hour[9] = 24*[1.]
    
    # SNAP 11 - Aviation
    tprof_mnth[10] = 12*[1.]
    tprof_week[10] =  7*[1.]
    tprof_hour[10] = 24*[1.]
    
    return tprof_hour, tprof_week, tprof_mnth


def get_local_time(utc_dt, timezone_name='Europe/Amsterdam'):
    """
    Convert UTC datetime to local datetime.
    """
    from zoneinfo import ZoneInfo
    from datetime import timezone

    if utc_dt.tzinfo is None:
        utc_dt = utc_dt.replace(tzinfo=timezone.utc)

    return utc_dt.astimezone(ZoneInfo(timezone_name))


def checkbounds(domainbounds, x, y):
    
    import numpy as np

    xminidx, xmaxidx, yminidx, ymaxidx = np.nan, np.nan, np.nan, np.nan

    # Check domain bounds
    do_proceed = True
    if domainbounds['xmin'] >= np.min(x):
        xminidx = np.argmax(x >= domainbounds['xmin'])
    else:
        do_proceed = False
        print("Requested lower bound x too low, request:", domainbounds['xmin'], "source file min:", np.min(x))

    if domainbounds['xmax'] <= np.max(x):
        xmaxidx = np.argmax(x >= domainbounds['xmax'])
    else:
        do_proceed = False
        print("Requested upper bound x too high, request:", domainbounds['xmax'], "source file max:", np.max(x))

    if domainbounds['ymin'] >= np.min(y):
        yminidx = np.argmax(y >= domainbounds['ymin'])
    else:
        do_proceed = False
        print("Requested lower bound y too low, request:", domainbounds['ymin'], "source file min:", np.min(y))

    if domainbounds['ymax'] <= np.max(y):
        ymaxidx = np.argmax(y >= domainbounds['ymax'])
    else:
        do_proceed = False
        print("Requested upper bound y too high, request:", domainbounds['ymax'], "source file max:", np.max(y))
    
    return do_proceed, int(xminidx), int(xmaxidx), int(yminidx), int(ymaxidx)




def writereademission(domainbounds, tstart, tend, tracers, sources, categos, show_log=False):
    
    import numpy as np
    import netCDF4 as netc
    import os
    import matplotlib.pyplot as plt
    from matplotlib import colors
    from datetime import datetime
    from datetime import timedelta
    
    """
    Read tracer sources of all categories, apply temporal profile and combine in 1 new file per tracer.
    
    Source files have to be in same directory and have to following form:
    <tracer>_<year>_<category>.nc
    In-file variables should be at least: 
    x-coordinates          (rank 1, size M) 
    y-coordinates          (rank 1, size N)
    tracer yearly emission (rank 2, size MxN) e.g. "co2", "ch4"
    
    """
    
    tstep = tstart
    dt = timedelta(hours=1)
    do_proceed = True
    hourly_emis = np.zeros([10, 10])  # dummy
    
    tprof_hour, tprof_week, tprof_mnth = loadsnap()
    
    while tstep <= tend and do_proceed:

        efs = tprof_mnth[:, tstep.month-1]*tprof_week[:, tstep.weekday()]*tprof_hour[:, tstep.hour]

        for itrac, tracer in enumerate(tracers):
            
            # Reset for new tracer
            hourly_emis[:] = 0
            
            tfname = "{}_emis_{}.nc".format(tracer, tstep.strftime("%Y%m%d%H%M"))
            if show_log:
                print("Target file: ", tfname)
            
            # -- Loop over source categories
            for isrc, sourcename in enumerate(sources[itrac]):
                
                icat = int(categos[itrac][isrc])

                # -- Read field
                sfname = "{}_{}_{}.nc".format(tracer, str(tstep.year), sourcename)
                if show_log:
                    print("Adding source: {:30s} (SNAP {:02d})".format(sfname, icat))

                if os.path.isfile(sfname):
                    sfobj  = netc.Dataset(sfname, 'r')
                    x = sfobj.variables['x'][:]
                    y = sfobj.variables['y'][:]
                    do_proceed, xminidx, xmaxidx, yminidx, ymaxidx = checkbounds(domainbounds, x, y)

                    if (itrac + isrc) == 0:
                        hourly_emis = np.zeros([xmaxidx-xminidx, ymaxidx-yminidx])

                    hourly_emis += sfobj.variables[tracer.upper()][xminidx:xmaxidx, yminidx:ymaxidx]*efs[icat-1]/(365*24)
                    sfobj.close()
                else:
                    do_proceed = False
                    print("File does not exist:", sfname)

            if do_proceed:
                # --- Create netCDF file
                tfobj = netc.Dataset(tfname, 'w')

                # --- Global attributes
                tfobj.title   = "100x100m " + tracer + " hourly emission map for " + tstep.strftime("%Y-%m-%d %H:%M") + "-" + (tstep+dt).strftime("%H:%M")
                tfobj.history = "Created: " + datetime.now().strftime("%d %b %Y")

                tfobj.description = "Total " + tracer + " emission from categories <list>"
                tfobj.valid       = "Valid from " + tstep.strftime("%Y-%m-%d %H:%M") + ' to ' + (tstep+dt).strftime("%Y-%m-%d %H:%M")
                tfobj.author      = 'A. Doyennel (VU)'
                tfobj.email       = 'a.doyennel@vu.nl'

                # -- Declaration of dimensions and variables
                dim_x = tfobj.createDimension('x', xmaxidx - xminidx)
                dim_y = tfobj.createDimension('y', ymaxidx - yminidx)

                var_x = tfobj.createVariable('x', 'f4', ('x', ))
                var_y = tfobj.createVariable('y', 'f4', ('y', ))
                var_e = tfobj.createVariable(tracer.lower(), 'f8', ('y', 'x', ))

                # -- Assigning values to variables
                var_x[:] = x[xminidx:xmaxidx]
                var_y[:] = y[yminidx:ymaxidx]
                var_e[:, :] = hourly_emis.T

                # -- Variable attributes
                var_x.units = 'Rijksdriehoekcoordinaat x in meters'
                var_y.units = 'Rijksdriehoekcoordinaat y in meters'
                var_e.units = 'Kilogram ' + tracer.lower() + ' per hour (kg hour-1)'

                tfobj.close()

                plt.figure(figsize=[12,8])
                plt.pcolormesh(x[xminidx:xmaxidx],y[yminidx:ymaxidx],hourly_emis.T, norm=colors.LogNorm(vmin=1e-1, vmax=1e4), cmap='Spectral_r')
                plt.text(915000,1045000,f'{tstep.strftime("%Y%m%d%H%M")}',fontsize=16)
                plt.colorbar(pad=.05, aspect=30)
                plt.tight_layout()
                plt.savefig("{}_emis_{}.png".format(tracer, tstep.strftime("%Y%m%d%H%M")))
                plt.close()

            tstep += dt




def writereademission_3d(inputncdir, output_dir, domainbounds, tstart, tend, tracers, sources, categos, tprof_hour, tprof_week, tprof_mnth, author, email, show_log=False):
    import numpy as np
    import netCDF4 as netc
    import os
    import matplotlib.pyplot as plt
    from matplotlib import colors
    from datetime import datetime
    from datetime import timedelta

    """
    Read tracer sources of all categories, apply temporal profile and combine in 1 new file per tracer.

    Source files have to be in same directory and have to following form:
    <tracer>_<year>_<category>.nc
    In-file variables should be at least: 
    x-coordinates          (rank 1, size X) 
    y-coordinates          (rank 1, size Y)
    z-coordinates          (rank 1, size Z)
    tracer yearly emission (rank 3, size X . Y . Z) e.g. "co2", "ch4"

    """

    tstep = tstart
    dt = timedelta(hours=1)
    do_proceed = True
    hourly_emis = np.zeros([10, 10])  # dummy

    #tprof_hour, tprof_week, tprof_mnth = loadsnap()

    while tstep <= tend and do_proceed:

        utc_t = tstep
        local_t = get_local_time(utc_t)

        efs = (
            tprof_mnth[:, local_t.month - 1]
            * tprof_week[:, local_t.weekday()]
            * tprof_hour[:, local_t.hour]
        )

        for itrac, tracer in enumerate(tracers):

            # Reset for new tracer
            hourly_emis[:] = 0

            tfname = "{}_emis_{}_3d.nc".format(tracer, utc_t.strftime("%Y%m%d%H%M"))
            if show_log:
                print("Target file: ", tfname)

            # -- Loop over source categories
            for isrc, sourcename in enumerate(sources[itrac]):

                icat = int(categos[itrac][isrc])

                # -- Read field
                sfname = "{}_{}_{}.nc".format(tracer, str(tstep.year), sourcename)
                if show_log:
                    print("Adding source: {:30s} (SNAP {:02d})".format(sfname, icat))

                if os.path.isfile(inputncdir + sfname):
                    sfobj = netc.Dataset(inputncdir + sfname, 'r')
                    rank = len(sfobj.variables[tracer.upper()].shape)
                    
                    if rank == 3:
                    
                        x = sfobj.variables['x'][:]
                        y = sfobj.variables['y'][:]
                        z = sfobj.variables['z'][:]

                        do_proceed, xminidx, xmaxidx, yminidx, ymaxidx, zminidx, zmaxidx = checkbounds_3d(domainbounds, x, y, z)

                        nx = xmaxidx - xminidx + 1
                        ny = ymaxidx - yminidx + 1
                        nz = zmaxidx - zminidx + 1
                        
                        if nx == len(x): 
                            sx = np.s_[:]
                        else:
                            sx = np.s_[xminidx:xmaxidx+1]
                            
                        if ny == len(y): 
                            sy = np.s_[:]
                        else:
                            sy = np.s_[yminidx:ymaxidx+1]         
                            
                        if nz == len(z): 
                            sz = np.s_[:]
                        else:
                            sz = np.s_[zminidx:zmaxidx+1]        
                            
                        if (itrac + isrc) == 0:
                            hourly_emis = np.zeros([nz, ny, nx])
                            
                        if show_log:
                            print(f"Field dimensions:        ({nz}, {ny}, {nx})")
                            print(f"Source file dimensions:  {sfobj.variables[tracer.upper()][sz, sy, sx].shape}")
                            print(f"Target array dimensions: {hourly_emis.shape}")

                        hourly_emis += sfobj.variables[tracer.upper()][sz, sy, sx] * efs[icat - 1] / (365 * 24)
                        
                    elif rank == 2:
                        
                        x = sfobj.variables['x'][:]
                        y = sfobj.variables['y'][:]

                        do_proceed, xminidx, xmaxidx, yminidx, ymaxidx = checkbounds(domainbounds, x, y)
                        
                        nx = xmaxidx - xminidx + 1
                        ny = ymaxidx - yminidx + 1
                        
                        if nx == len(x): 
                            sx = np.s_[:]
                        else:
                            sx = np.s_[xminidx:xmaxidx+1]
                            
                        if ny == len(y): 
                            sy = np.s_[:]
                        else:
                            sy = np.s_[yminidx:ymaxidx+1]
  
                        if (itrac + isrc) == 0:
                            print("ERROR: Initializing with 2D! Reconsider code!")
                            break

                        hourly_emis[0] += sfobj.variables[tracer.upper()][sx, sy].T * efs[icat - 1] / (365 * 24)    
                        
                    sfobj.close()
                else:
                    do_proceed = False
                    print("File does not exist:", sfname)

            if do_proceed:
                # --- Create netCDF file
                tfobj = netc.Dataset(output_dir + tfname, 'w')


                t0 = utc_t.strftime("%Y-%m-%d %H:%M")
                t1 = (utc_t + dt).strftime("%H:%M")

                
                dx = float(x[1]-x[0])
                dy = float(y[1]-y[0])
                
                # --- Global attributes
                tfobj.title = f"{dx}x{dy}m {tracer.upper()} hourly emissions map for {t0}-{t1}"
                tfobj.history = "Created: " + datetime.now().strftime("%d %b %Y")

                tfobj.description = f"Total {tracer.upper()} emission from categories {sources}"
                tfobj.valid = ("Valid from " + utc_t.strftime("%Y-%m-%d %H:%M") + " to " + (utc_t + dt).strftime("%Y-%m-%d %H:%M"))
                tfobj.author = author
                tfobj.email = email

                # -- Declaration of dimensions and variables
                dim_x = tfobj.createDimension('x', nx)
                dim_y = tfobj.createDimension('y', ny)
                dim_z = tfobj.createDimension('z', nz)

                var_x = tfobj.createVariable('x', 'f4', ('x',))
                var_y = tfobj.createVariable('y', 'f4', ('y',))
                var_z = tfobj.createVariable('z', 'f4', ('z',))
                var_e = tfobj.createVariable(tracer.lower(), 'f8', ('z', 'y', 'x',))

                # -- Assigning values to variables
                var_x[:] = x[sx]
                var_y[:] = y[sy]
                var_z[:] = z[sz]
                var_e[:, :, :] = hourly_emis

                # -- Variable attributes
                var_x.units = 'Rijksdriehoekcoordinaat x in meters'
                var_y.units = 'Rijksdriehoekcoordinaat y in meters'
                var_z.units = 'Hoogte z in meters (middelpunt gridbox)'
                var_e.units = 'kg hour-1'
                #var_e.units = 'Kilogram ' + tracer.lower() + ' per hour (kg hour-1)'

                tfobj.close()

                # plt.figure(figsize=[12, 8])
                # plt.pcolormesh(x[xminidx:xmaxidx], y[yminidx:ymaxidx], hourly_emis[:,:,5].T,
                #                norm=colors.LogNorm(vmin=1e-1, vmax=1e4), cmap='Spectral_r')
                # plt.text(915000, 1045000, f'{tstep.strftime("%Y%m%d%H%M")}', fontsize=16)
                # plt.colorbar(pad=.05, aspect=30)
                # plt.tight_layout()
                # plt.savefig("{}_emis_{}.png".format(tracer, tstep.strftime("%Y%m%d%H%M")))
                # plt.close()

            tstep += dt


def extract_multicat(sourcefile, subcats=[], cats=[], year=-1, xmin=-1, xmax=-1, dx=-1, ymin=-1, ymax=-1, dy=-1):

    import pandas as pd
    import numpy as np

    """
    Extract and combine specified emission category and subcategories from ER 1x1km file.
    """

    if any(t == -1 for t in [year, xmin, xmax, dx, ymin, ymax, dy]):
        print('Missing input found')
        print('jaar {} xmin {} xmax {} dx {} ymin {} ymax {} dy {} '.format(year, xmin, xmax, dx, ymin, ymax, dy))
        return

    if (len(cats) + len(subcats)) == 0:
        print('No categories or subcategories found')
        return

    # ===== Prepare spatial data ====================
    xmin = int(xmin)
    xmax = int(xmax)
    dx = int(dx)
    ymin = int(ymin)
    ymax = int(ymax)
    dy = int(dy)

    x = range(xmin, xmax + 1, dx)
    y = range(ymin, ymax + 1, dy)

    # Construct numpy array based on required extent
    data = np.zeros([(xmax - xmin) // dx, (ymax - ymin) // dy])

    # =============================================
    df_main = pd.read_csv(sourcefile, delimiter=';')  # Read in file

    if len(cats) > 0:
        for cat in cats:
            df = df_main.loc[df_main['DOELGROEP'] == cat]  # Select doelgroep
            data = df2numpy(df, data, year, xmin, xmax, dx, ymin, ymax, dy)

    if len(subcats) > 0:
        for subcat in subcats:
            df = df_main.loc[df_main['SUBDOELGROEP'] == subcat]  # Select subdoelgroep
            data = df2numpy(df, data, year, xmin, xmax, dx, ymin, ymax, dy)

    return data, np.array(x), np.array(y)




def shrink(data, rows, cols):
    return data.reshape(rows, data.shape[0]//rows, cols, data.shape[1]//cols).sum(axis=1).sum(axis=2)


def tile_array(a, b0, b1):
    
    from numpy.lib.stride_tricks import as_strided
    
    r, c = a.shape                                    # number of rows/columns
    rs, cs = a.strides                                # row/column strides 
    x = as_strided(a, (r, b0, c, b1), (rs, 0, cs, 0)) # view a as larger 4D array
    return x.reshape(r*b0, c*b1)                      # create new 2D array



def load_ERmap(df, field, x, y, year):
    import numpy as np
    
    year = str(year)
    
    df_5km = df[df['CODE_GEBIED'].str.contains('km5')]
    df_1km = df[df['CODE_GEBIED'].str.contains('km1')]
    
    def convert_to_float(value):
        if isinstance(value, str):
            value = value.replace(',', '.')
        return float(value)
    
    if len(df_5km) > 0:
        for index, row in df_5km.iterrows():
            mapemission = convert_to_float(row['EMISSIE_KG'])
            
            xidx = np.argmax(x >= 1e3 * float(row['CODE_GEBIED'][1:4]))
            yidx = np.argmax(y >= 1e3 * float(row['CODE_GEBIED'][5:8]))
            field[xidx:xidx + 5, yidx:yidx + 5] += mapemission / 25
    
    if len(df_1km) > 0:
        for index, row in df_1km.iterrows():
            mapemission = convert_to_float(row['EMISSIE_KG'])
    
            xidx = np.argmax(x >= 1e3 * float(row['CODE_GEBIED'][0:3]))
            yidx = np.argmax(y >= 1e3 * float(row['CODE_GEBIED'][3:6]))
            field[xidx, yidx] += mapemission

import numpy as np

def add_subtract(point_emis, mode):
    return point_emis if mode == 'add' else  point_emis*0 #-point_emis #(if point sources are not included to area emissions, use point_emis*0 !!!!)

def join_map_points(df, field, x=[], y=[], z=[], mode=None, dim=2, xloc='XCOORD', yloc='YCOORD', zloc='HOOGTE', emisname='EMISSIE'):

    if mode not in ['add', 'subtract']:
        raise ValueError('Invalid mode, use "add" or "subtract"')
    if dim not in [2, 3]:
        raise ValueError('Invalid dimension, use 2 or 3')

    # -------------------------
    # CLEAN INPUT
    # -------------------------
    df = df.dropna(subset=[xloc, yloc, emisname]).copy()

    n_in = len(df)
    n_mapped = 0
    n_skipped = 0
    grid_updates = 0

    # -------------------------
    # MAIN LOOP
    # -------------------------
    for _, row in df.iterrows():

        reason_skipped = None

        x_point = float(row[xloc])
        y_point = float(row[yloc])
        emis = float(row[emisname])

        # --- emission validity ---
        if not np.isfinite(emis):
            n_skipped += 1
            continue

        # --- domain check ---
        if not (x[0] <= x_point <= x[-1] and y[0] <= y_point <= y[-1]):
            n_skipped += 1
            continue

        # --- indices ---
        xidx = np.searchsorted(x, x_point) - 1
        yidx = np.searchsorted(y, y_point) - 1

        # --- clamp safety (IMPORTANT FIX) ---
        if xidx < 0 or xidx >= len(x) - 1 or yidx < 0 or yidx >= len(y) - 1:
            n_skipped += 1
            continue

        # -------------------------
        # VERTICAL DIMENSION
        # -------------------------
        if dim == 3:
            if zloc not in row or not np.isfinite(row[zloc]):
                print('WARNING: point source file has no stack height info -> point sources are added to the lowest z')
                z_point = z[0]
            else:
                z_point = float(row[zloc])

            zidx = np.searchsorted(z, z_point) - 1
            zidx = max(min(zidx, len(z) - 2), 0)
        else:
            zidx = None

        # -------------------------
        # POINT SOURCE
        # -------------------------
        if row['TYPE'] == 'P':

            if dim == 2:
                field[xidx, yidx] += add_subtract(emis, mode)
            else:
                field[xidx, yidx, zidx] += add_subtract(emis, mode)

            grid_updates += 1
            n_mapped += 1

        # -------------------------
        # AREA SOURCE
        # -------------------------
        else:

            length = row['LENGTE']
            width = row['BREEDTE']

            x_upper, x_lower = x_point + length / 2, x_point - length / 2
            y_upper, y_lower = y_point + width / 2, y_point - width / 2

            # skip if outside domain
            if x_upper < x[0] or x_lower > x[-1] or y_upper < y[0] or y_lower > y[-1]:
                n_skipped += 1
                continue

            xidx_0 = np.searchsorted(x, x_lower) - 1
            xidx_1 = np.searchsorted(x, x_upper) - 1
            yidx_0 = np.searchsorted(y, y_lower) - 1
            yidx_1 = np.searchsorted(y, y_upper) - 1

            point_used = False

            for ix in range(xidx_0, xidx_1 + 1):
                for iy in range(yidx_0, yidx_1 + 1):

                    fracx = (min(x_upper, x[ix + 1]) - max(x[ix], x_lower)) / length
                    fracy = (min(y_upper, y[iy + 1]) - max(y[iy], y_lower)) / width

                    if fracx > 0 and fracy > 0:

                        if dim == 2:
                            field[ix, iy] += add_subtract(emis * fracx * fracy, mode)
                        else:
                            field[ix, iy, zidx] += add_subtract(emis * fracx * fracy, mode)

                        grid_updates += 1
                        point_used = True

            if point_used:
                n_mapped += 1
            else:
                n_skipped += 1

    # -------------------------
    # FINAL LOG (CORRECT)
    # -------------------------
    print("--------------------------------------------------")
    print(f"INPUT POINTS: {n_in}")
    print(f"MAPPED POINTS (P + O): {n_mapped}")
    print(f"SKIPPED POINTS: {n_skipped}")
    print(f"GRID CELL UPDATES (not points): {grid_updates}")
    print("--------------------------------------------------")
                        

def create3Dfromscratch(x, y, z, zh, projname, snap, pointfile, verbose=False, savepoints=False):

    total_emis = 0

    import numpy as np
    import pandas as pd

    xmax = np.max(x)
    xmin = np.min(x)
    dx = x[1] - x[0]

    ymax = np.max(y)
    ymin = np.min(y)
    dy = y[1] - y[0]

    emis_points = np.zeros([len(x), len(y), len(z)])

    df = pd.read_csv(pointfile, delimiter=',', decimal=".", encoding="utf-8")

    xlocname = f'EMISSION_LOCATION_X_{projname}'
    ylocname = f'EMISSION_LOCATION_Y_{projname}'
    zlocname = 'EMISSION_LOCATION_Z'

    xextname = 'EMISSION_LOCATION_X_RANGE'
    yextname = 'EMISSION_LOCATION_Y_RANGE'
    zextname = 'EMISSION_LOCATION_Z_RANGE'

    df = df.loc[(df[xlocname] < xmax + dx) &
                (df[xlocname] >= xmin) &
                (df[ylocname] < ymax + dy) &
                (df[ylocname] >= ymin) &
                (df['SNAP'] == snap)]

    print(snap, np.shape(df))
    iirow = 0
    for irow, row in df.iterrows():

        iirow += 1
        point_x = row[xlocname]
        point_y = row[ylocname]
        point_z = row[zlocname]

        point_x_ext = row[xextname]
        point_y_ext = row[yextname]
        point_z_ext = row[zextname]

        point_e = row['EMISSIE']

        point_x_idx_low  = np.argmax(x >  (point_x - point_x_ext / 2)) - 1
        point_x_idx_high = np.argmax(x >= (point_x + point_x_ext / 2)) - 1

        point_y_idx_low  = np.argmax(y >  (point_y - point_y_ext / 2)) - 1
        point_y_idx_high = np.argmax(y >= (point_y + point_y_ext / 2)) - 1

        point_z_idx_low  = np.argmax(zh >  (point_z - point_z_ext / 2)) - 1
        point_z_idx_high = np.argmax(zh >= (point_z + point_z_ext / 2)) - 1

        xn = max(point_x_idx_high - point_x_idx_low + 1, 1)
        yn = max(point_y_idx_high - point_y_idx_low + 1, 1)
        zn = max(point_z_idx_high - point_z_idx_low + 1, 1)

        emis_points[point_x_idx_low: point_x_idx_low + xn,
                    point_y_idx_low: point_y_idx_low + yn,
                    point_z_idx_low: point_z_idx_low + zn] += point_e / (xn * yn * zn)

        # print(xn,yn,zn,xn*yn*zn, np.sum(np.shape(emis_points[point_x_idx_low: point_x_idx_low + xn,
        #                                  point_y_idx_low: point_y_idx_low + yn,
        #                                  point_z_idx_low: point_z_idx_low + zn])))
        total_emis += point_e

        if verbose: print(iirow, row['EMISSIEOORZAAK'], point_x, point_e, xn, yn, zn)
    if savepoints: 
        df.to_csv(f'co2_2017_SNAP_{snap}_points.csv')
        print(f'create3Dfromscratch: Point source emissions written to file: co2_2017_SNAP_{snap}_points.csv')
            
    print(total_emis, np.sum(emis_points))
    return emis_points
    
    
    
def pointsto3D(df, xt, yt, zt, proj='RD', verbose=False):
    
    import numpy as np
    import pandas as pd
    import pyproj
    
    """
    Input:
    df    - Dataframe with point emissions
    xt,yt,zt - (Numpy) arrays with gridcell edges, nb: length is number of gridcells + 1
    
    Output:
    
    """
    
    emis = np.zeros([len(xt)-1,len(yt)-1,len(zt)-1])
        
    for irow, row in df.iterrows():
        
        # EMISSIE
        point_e = row['EMISSIE']
        
        # POINT SOURCE COORDINATES (AND EXTENSION)
        if proj == 'HARM':
            point_x = row['XCO_EMISSIEPUNT_HARM'] 
            point_y = row['YCO_EMISSIEPUNT_HARM']
        else:
            point_x = row['XCO_EMISSIEPUNT'] 
            point_y = row['YCO_EMISSIEPUNT']
        
        if not point_x and point_y:
            if proj == 'HARM':
                point_x = row['XCO_BEDRIJF_HARM'] 
                point_y = row['YCO_BEDRIJF_HARM']
            else:
                point_x = row['XCO_BEDRIJF'] 
                point_y = row['YCO_BEDRIJF']
        
        point_x_ext = np.nan_to_num(row['LENGTE'] )
        point_y_ext = np.nan_to_num(row['BREEDTE'])
        
        # POINT SOURCE HEIGHT
        point_z = np.nan_to_num(row['HOOGTE'])
            
        if 'Vliegverkeer' in row['EMISSIEOORZAAK']:
            if point_z == 0:
                point_z_ext = 0
                
            elif point_z < 200:
                point_z_ext = 50
                
            elif point_z < 500:
                point_z_ext = 100
                
            elif point_z < 1000:
                point_z_ext = 250
                
            else:
                point_z_ext = 0
                print("Vliegverkeer altitude >1000", point_z, row['EMISSIEOORZAAK'])
        else:
            point_z_ext = 0
        
        point_x_idx_low  = np.argmax(xt> (point_x - point_x_ext/2))-1
        point_x_idx_high = np.argmax(xt>=(point_x + point_x_ext/2))-1
        xn = point_x_idx_high - point_x_idx_low + 1
        
        point_y_idx_low  = np.argmax(yt> (point_y - point_y_ext/2))-1
        point_y_idx_high = np.argmax(yt>=(point_y + point_y_ext/2))-1
        yn = point_y_idx_high - point_y_idx_low + 1
            
        point_z_idx_low  = np.argmax(zt> (point_z - point_z_ext/2))-1
        point_z_idx_high = np.argmax(zt>=(point_z + point_z_ext/2))-1
        zn = point_z_idx_high - point_z_idx_low + 1
        
        emis[point_x_idx_low : point_x_idx_low + xn, 
             point_y_idx_low : point_y_idx_low + yn, 
             point_z_idx_low : point_z_idx_low + zn] += point_e/(xn*yn*zn)
        
    return emis    
    


def checkbounds_3d(domainbounds, x, y, z):
    
    import numpy as np
    import xarray as xr

    xminidx, xmaxidx, yminidx, ymaxidx, zminidx, zmaxidx = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    # Check domain bounds
    # --- XMIN ------------------------------------------------------------------------------------------------
    do_proceed = True
    if domainbounds['xmin'] >= np.min(x):
        xminidx = np.argmax(x >= domainbounds['xmin'])
    else:
        do_proceed = False
        print("Requested lower bound x too low, request:", domainbounds['xmin'], "source file min:", np.min(x))
    
    # --- XMAX ------------------------------------------------------------------------------------------------
    if domainbounds['xmax'] <= np.max(x):
        xmaxidx = np.argmax(x >= domainbounds['xmax'])
    else:
        do_proceed = False
        print("Requested upper bound x too high, request:", domainbounds['xmax'], "source file max:", np.max(x))
    
    # --- YMIN ------------------------------------------------------------------------------------------------
    if domainbounds['ymin'] >= np.min(y):
        yminidx = np.argmax(y >= domainbounds['ymin'])
    else:
        do_proceed = False
        print("Requested lower bound y too low, request:", domainbounds['ymin'], "source file min:", np.min(y))

    # --- YMAX ------------------------------------------------------------------------------------------------
    if domainbounds['ymax'] <= np.max(y):
        ymaxidx = np.argmax(y >= domainbounds['ymax'])
    else:
        do_proceed = False
        print("Requested upper bound y too high, request:", domainbounds['ymax'], "source file max:", np.max(y))
    
    # --- ZMIN ------------------------------------------------------------------------------------------------
    # Since this is about altitude, it does not make sense to check minimum, since we always work from the ground up.
    zminidx = np.argmax(z >= domainbounds['zmin'])

    
    # --- ZMAX ------------------------------------------------------------------------------------------------
    if domainbounds['zmax'] <= np.max(z):
        zmaxidx = np.argmax(z >= domainbounds['zmax'])
    else:
        do_proceed = False
        print("Requested upper bound z too high, request:", domainbounds['zmax'], "source file max:", np.max(z))
    
    print( xminidx, xmaxidx, yminidx, ymaxidx, zminidx, zmaxidx )   
    return do_proceed, int(xminidx), int(xmaxidx), int(yminidx), int(ymaxidx), int(zminidx), int(zmaxidx)


  

def df2numpy(df, data, year, xmin, xmax, dx, ymin, ymax, dy):

    import numpy as np

    for idx, content in df.iterrows():

        location = content['CODE_GEBIED']
        emis = content['EMISSIE_KG']
        emis = float(emis.replace(',','.'))
        
        
        emis = content[str(year)]

        if location[-3:] == 'km1':
            xidx = int((int(location[0:3]) * 1e3 - xmin) // dx)
            yidx = int((int(location[3:6]) * 1e3 - ymin) // dy)

            if (xidx >= 0) and (yidx >= 0) and (xidx < data.shape[0]) and (yidx < data.shape[1]):
                data[xidx, yidx] += emis

    return data
    
    
    
    
    
def downscale(sourcefield, source_dx, target_dx, pos_only=True):
    import numpy as np

    if pos_only:
        sourcefield = np.array(sourcefield)
        sourcefield[sourcefield == -99997] = 0

    def rebin(a, shape):
        sh = shape[0], a.shape[0] // shape[0], shape[1], a.shape[1] // shape[1]
        return a.reshape(sh).sum(-1).sum(1)

    if np.mod(target_dx, source_dx) != 0:
        print(
            'ERROR: Target resolution (dx={}) not a multiple of source resolution (dx={})'.format(target_dx, source_dx))
        return

    targetshape = (int(np.ceil(np.shape(sourcefield)[0] * source_dx / target_dx)),
                   int(np.ceil(np.shape(sourcefield)[1] * source_dx / target_dx)))
    targetfield = rebin(sourcefield, targetshape)

    return targetfield
    
    
    
    
def lowonhigh(hires, lores, do_fill=False, fillvalue=999.):
    import numpy as np

    shapelores = np.shape(lores)
    shapehires = np.shape(hires)
    ratio = np.array(shapehires) / np.array(shapelores)

    if np.all(np.mod(ratio, 1) == 0):
#         print('Ratio ok!', ratio)
        ratio = ratio.astype(int)
    else:
        print('Ratio not a whole number!', ratio)
        return

    target = np.zeros(shapehires)

    for ix in range(shapelores[0]):
        for iy in range(shapelores[1]):
            value = lores[ix, iy]
            if do_fill and (value <= 0):
                value = fillvalue
            target[ix * ratio[0]:(ix + 1) * ratio[0], iy * ratio[1]:(iy + 1) * ratio[1]] = value

    return target
    
    
    

def extractpoints(sourcefile, fldname, xmin, xmax, dx, ymin, ymax, dy, zmin, zmax, dz, pos_only=True, verbose=False):

    import pandas as pd
    import numpy as np

    """

    """

    xmin, xmax, dx = int(xmin), int(xmax), int(dx)
    ymin, ymax, dy = int(ymin), int(ymax), int(dy)
    zmin, zmax, dz = int(zmin), int(zmax), int(dz)

    x = range(xmin, xmax, dx)
    y = range(ymin, ymax, dy)
    z = range(zmin, zmax, dz)

    df = pd.read_csv(sourcefile, delimiter=';', decimal=",")  # Read in file

    # Construct numpy array based on required extent
    fillval = 0
    data = np.ones([(xmax - xmin) // dx, (ymax - ymin) // dy, (zmax - zmin) // dz]) * fillval

    for idx, content in df.iterrows():

        xloc = content['XCO_EMISSIEPUNT']
        yloc = content['YCO_EMISSIEPUNT']
        zloc = content['HOOGTE']
        emis = content[fldname]
        label = content['TYPE']

        if np.isnan(xloc) or np.isnan(yloc) or np.isnan(zloc) or np.isnan(emis) or not label == 'P':
            if verbose: print(label, xloc, yloc, zloc, emis)

        else:

            xloc = int(xloc)
            yloc = int(yloc)
            zloc = int(zloc)
            emis = float(emis)

            xidx = (xloc - xmin) // dx
            yidx = (yloc - ymin) // dy
            zidx = (zloc - zmin) // dz

            if (xidx >= 0) and (yidx >= 0) and (zidx >= 0) and (xidx < len(x)) and (yidx < len(y) and (zidx < len(z))):
                if pos_only:
                    data[xidx, yidx, zidx] = data[xidx, yidx, zidx] + max(0, emis)
                else:
                    data[xidx, yidx, zidx] = data[xidx, yidx, zidx] + emis
            else:
                print('Outside domain bounds:', 'xloc', xidx, 'yloc', yloc, 'zloc', zloc, 'data', emis)

    return data, np.array(x), np.array(y), np.array(z)



def calcppm(emis_input, rho_air=1.17e3, dx=100, dy=100, dz=40):
    m_co2 = 44
    m_air = 28.94

    emis_year = emis_input * 1e3  # grams   year-1
    emis_min = emis_year / (365 * 24 * 60)  # grams minute-1
    emis_mol = emis_min / m_co2  # mol   minute-1

    dv = dx * dy * dz

    g_air = rho_air * dv  # g gridbox-1
    air_mol = g_air / m_air  # mol gridbox-1

    emis_ppm = emis_mol / (emis_mol + air_mol) * 1e6

    return emis_ppm
    
    
    



def extractpoints2(df, xmin, xmax, dx, ymin, ymax, dy, zmin, zmax, dz, pos_only=True, verbose=False):
    
    import numpy as np

    """

    """
    xmin, xmax, dx = int(xmin), int(xmax), int(dx)
    ymin, ymax, dy = int(ymin), int(ymax), int(dy)
    zmin, zmax, dz = int(zmin), int(zmax), int(dz)

    x = range(xmin, xmax, dx)
    y = range(ymin, ymax, dy)
    z = range(zmin, zmax, dz)

    # Construct numpy array based on required extent
    fillval = 0
    data = np.ones([(xmax - xmin) // dx, (ymax - ymin) // dy, (zmax - zmin) // dz]) * fillval

    nfills = 0
    for idx, content in df.iterrows():

        xloc = content['XCO_EMISSIEPUNT']
        yloc = content['YCO_EMISSIEPUNT']
        zloc = content['HOOGTE']
        emis = content['EMISSIE']
        label = content['TYPE']

        if ('RWZI' in content['EMISSIEOORZAAK'] or 'SBI 37' in content['EMISSIEOORZAAK']) :
            xloc = content['XCO_BEDRIJF']
            yloc = content['YCO_BEDRIJF']
            zloc = 0

        if (np.isnan(xloc) or np.isnan(yloc) or np.isnan(zloc) or np.isnan(emis)):

            if verbose: print(content['EMISSIEOORZAAK'], label, xloc, yloc, zloc, emis)

        else:

            xloc = int(xloc)
            yloc = int(yloc)
            zloc = int(zloc)
            emis = float(emis)

            xidx = (xloc - xmin) // dx
            yidx = (yloc - ymin) // dy
            zidx = (zloc - zmin) // dz

            if (xidx >= 0) and (yidx >= 0) and (zidx >= 0) and (xidx < len(x)) and (yidx < len(y) and (zidx < len(z))):
                nfills = nfills+1
                if pos_only:
                    data[xidx, yidx, zidx] = data[xidx, yidx, zidx] + max(0, emis)
                else:
                    data[xidx, yidx, zidx] = data[xidx, yidx, zidx] + emis
            else:
                if verbose: print('Outside domain bounds:', 'xloc', xidx, 'yloc', yloc, 'zloc', zloc, 'data', emis)
    print(nfills)
    return data, np.array(x), np.array(y), np.array(z)




def extractpoints2_withtransform(df, xmin, xmax, dx, ymin, ymax, dy, zmin, zmax, dz, pos_only=True, verbose=False):
    import numpy as np
    import pyproj as pyproj
    """

    """
    xmin, xmax, dx = int(xmin), int(xmax), int(dx)
    ymin, ymax, dy = int(ymin), int(ymax), int(dy)
    zmin, zmax, dz = int(zmin), int(zmax), int(dz)

    x = range(xmin, xmax, dx)
    y = range(ymin, ymax, dy)
    z = range(zmin, zmax, dz)

    # Construct numpy array based on required extent
    fillval = 0
    data = np.ones([(xmax - xmin) // dx, (ymax - ymin) // dy, (zmax - zmin) // dz]) * fillval

    nfills = 0
    for idx, content in df.iterrows():

        xloc = content['XCO_EMISSIEPUNT']
        yloc = content['YCO_EMISSIEPUNT']
        zloc = content['HOOGTE']
        emis = content['EMISSIE']
        label = content['TYPE']

        if ('RWZI' in content['EMISSIEOORZAAK'] or 'SBI 37' in content['EMISSIEOORZAAK']):
            xloc = content['XCO_BEDRIJF']
            yloc = content['YCO_BEDRIJF']
            zloc = 0

        if (np.isnan(xloc) or np.isnan(yloc) or np.isnan(zloc) or np.isnan(emis)):

            if verbose: print(content['EMISSIEOORZAAK'], label, xloc, yloc, zloc, emis)

        else:

            proj_RD   = pyproj.Proj('+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.417,50.3319,465.552,-0.398957,0.343988,-1.8774,4.0725 +units=m +no_defs')
            proj_HARM = pyproj.Proj('+proj=lcc +lat_1=52.500000 +lat_2=52.500000 +lat_0=52.500000 +lon_0=.000000 +k_0=1.0 +x_0=649536.512574 +y_0=1032883.739533 +a=6371220.000000 +b=6371220.000000')

            xloc, yloc = pyproj.transform(proj_RD, proj_HARM, xloc, yloc)

            xloc = int(xloc)
            yloc = int(yloc)
            zloc = int(zloc)
            emis = float(emis)

            xidx = (xloc - xmin) // dx
            yidx = (yloc - ymin) // dy
            zidx = (zloc - zmin) // dz

            if (xidx >= 0) and (yidx >= 0) and (zidx >= 0) and (xidx < len(x)) and (yidx < len(y) and (zidx < len(z))):
                nfills = nfills + 1
                if pos_only:
                    data[xidx, yidx, zidx] = data[xidx, yidx, zidx] + max(0, emis)
                else:
                    data[xidx, yidx, zidx] = data[xidx, yidx, zidx] + emis
            else:
                if verbose: print('Outside domain bounds:', 'xloc', xloc, 'yloc', yloc, 'zloc', zloc, 'data', emis)
    print(nfills)
    return data, np.array(x), np.array(y), np.array(z)




def create3Dfrom2D(z, projname, snap, emissionfile, pointfile):

    import numpy as np
    import netCDF4 as netc
    import pandas as pd

    fobj = netc.Dataset(emissionfile, 'r')
    x = fobj.variables['x'][:]
    y = fobj.variables['y'][:]
    emis_orig = fobj.variables['CO2'][:]
    fobj.close()

    xmax = np.max(x)
    xmin = np.min(x)
    dx = x[1] - x[0]

    ymax = np.max(y)
    ymin = np.min(y)
    dy = y[1] - y[0]

    emis_points = np.zeros([len(x), len(y), len(z)])

    df = pd.read_csv(pointfile, delimiter=',', decimal=".", encoding="utf-8")

    xlocname = f'EMISSION_LOCATION_X_{projname}'
    ylocname = f'EMISSION_LOCATION_Y_{projname}'
    zlocname = 'EMISSION_LOCATION_Z'

    xextname = 'EMISSION_LOCATION_X_RANGE'
    yextname = 'EMISSION_LOCATION_Y_RANGE'
    zextname = 'EMISSION_LOCATION_Z_RANGE'

    df = df.loc[(df[xlocname] < xmax + dx) &
                (df[xlocname] >= xmin) &
                (df[ylocname] < ymax + dy) &
                (df[ylocname] >= ymin) &
                (df['SNAP'] == snap)]

    print(np.shape(df))
    iirow = 0
    for irow, row in df.iterrows():
        iirow += 1
        point_x = row[xlocname]
        point_y = row[ylocname]
        point_z = row[zlocname]

        point_x_ext = row[xextname]
        point_y_ext = row[yextname]
        point_z_ext = row[zextname]

        point_e = row['EMISSIE']

        point_x_idx_low = np.argmax(x > (point_x - point_x_ext / 2)) - 1
        point_x_idx_high = np.argmax(x >= (point_x + point_x_ext / 2)) - 1
        xn = point_x_idx_high - point_x_idx_low + 1

        point_y_idx_low = np.argmax(y > (point_y - point_y_ext / 2)) - 1
        point_y_idx_high = np.argmax(y >= (point_y + point_y_ext / 2)) - 1
        yn = point_y_idx_high - point_y_idx_low + 1

        point_z_idx_low = np.argmax(z > (point_z - point_z_ext / 2)) - 1
        point_z_idx_high = np.argmax(z >= (point_z + point_z_ext / 2)) - 1
        zn = max(point_z_idx_high - point_z_idx_low + 1, 1)

        emis_points[point_x_idx_low: point_x_idx_low + xn,
        point_y_idx_low: point_y_idx_low + yn,
        point_z_idx_low: point_z_idx_low + zn] += point_e / (xn * yn * zn)

        print(iirow, row['EMISSIEOORZAAK'], point_x, point_e, xn, yn, zn)

    return x, y, emis_orig, emis_points




def regrescat(dataframe, DoVerbose = False, v_order = 1, t_order = 0, h_order = 1, h_shape = 'log'):
    from sklearn.preprocessing import PolynomialFeatures
    import statsmodels.api as sm
    
    x_all = np.log10(dataframe.EMISSIE.     values)
    new_x = np.linspace(np.min(x_all),np.max(x_all),num=20)
    
    # VOLUMESTROOM
    dataframe_v = dataframe.loc[(dataframe['VOLUMESTROOM'] > 0.)]
    x = np.log10(dataframe_v.EMISSIE.     values)
    v = np.log10(dataframe_v.VOLUMESTROOM.values)
        
    x = PolynomialFeatures(v_order).fit_transform(x.reshape(-1,1))
        
    reg_v  = sm.OLS(v,x).fit()
    poly_v = np.poly1d(np.flipud(reg_v.params))
    new_v  = poly_v(new_x) 
    
    # TEMPERATUUR
    dataframe_t = dataframe.loc[(dataframe['TEMPERATUUR'] > 0.)]
    x = np.log10(dataframe_t.EMISSIE.     values)
    t =          dataframe_t.TEMPERATUUR. values
    
    x = PolynomialFeatures(t_order).fit_transform(x.reshape(-1,1))
    
    reg_t  = sm.OLS(t,x).fit()
    poly_t = np.poly1d(np.flipud(reg_t.params))
    new_t  = poly_t(new_x)      
    
    # HOOGTE
    dataframe_h = dataframe.loc[(dataframe['HOOGTE'] > 0.)]
    x = np.log10(dataframe_h.EMISSIE.     values)
    if h_shape == 'log':
        h = np.log10(dataframe_h.HOOGTE.     values)
    elif h_shape == 'exp':
        h = np.log(dataframe_h.HOOGTE.     values)
    else:
        h =          dataframe_h.HOOGTE.      values

    x = PolynomialFeatures(h_order).fit_transform(x.reshape(-1,1))
    
    if DoVerbose:
        print(x)
        
    reg_h  = sm.OLS(h,x,hasconst=True).fit()
    poly_h = np.poly1d(np.flipud(reg_h.params))
    new_h  = poly_h(new_x) 
    
    if h_shape == 'log':
        new_h = 10**new_h
    elif h_shape == 'exp': 
        new_h = np.exp(new_h)

    return 10**new_x, 10**new_v, new_t, new_h



def gapfill(dataframe, DoVerbose=False,
            v_order=1, t_order=0,
            h_order=1, h_shape='log',
            a_order=1, area_shape='log'):

    from sklearn.preprocessing import PolynomialFeatures
    import statsmodels.api as sm
    import numpy as np

    # ================================================================
    # 1. VOLUMESTROOM
    # ================================================================
    dataframe_v  = dataframe.loc[dataframe['VOLUMESTROOM'] > 0]
    dataframe_v0 = dataframe.loc[dataframe['VOLUMESTROOM'] == 0]

    x = np.log10(dataframe_v.EMISSIE.values)
    v = np.log10(dataframe_v.VOLUMESTROOM.values)

    x_poly = PolynomialFeatures(v_order).fit_transform(x.reshape(-1, 1))
    reg_v  = sm.OLS(v, x_poly).fit()
    poly_v = np.poly1d(np.flipud(reg_v.params))

    replace_v = []
    for idx, row in dataframe_v0.iterrows():
        replace_v.append([idx, 10**poly_v(np.log10(row['EMISSIE']))])

    del dataframe_v, dataframe_v0

    # ================================================================
    # 2. TEMPERATUUR
    # ================================================================
    dataframe_t  = dataframe.loc[dataframe['TEMPERATUUR'] > 0]
    dataframe_t0 = dataframe.loc[dataframe['TEMPERATUUR'] == 0]

    x = np.log10(dataframe_t.EMISSIE.values)
    t = dataframe_t.TEMPERATUUR.values

    x_poly = PolynomialFeatures(t_order).fit_transform(x.reshape(-1, 1))
    reg_t  = sm.OLS(t, x_poly).fit()
    poly_t = np.poly1d(np.flipud(reg_t.params))

    replace_t = []
    for idx, row in dataframe_t0.iterrows():
        replace_t.append([idx, poly_t(np.log10(row['EMISSIE']))])

    del dataframe_t, dataframe_t0

    # ================================================================
    # 3. HOOGTE
    # ================================================================
    dataframe_h  = dataframe.loc[dataframe['HOOGTE'] > 0]
    dataframe_h0 = dataframe.loc[dataframe['HOOGTE'] == 0]

    x = np.log10(dataframe_h.EMISSIE.values)

    if h_shape == 'log':
        h = np.log10(dataframe_h.HOOGTE.values)
    elif h_shape == 'exp':
        h = np.log(dataframe_h.HOOGTE.values)
    else:
        h = dataframe_h.HOOGTE.values

    x_poly = PolynomialFeatures(h_order).fit_transform(x.reshape(-1, 1))
    reg_h  = sm.OLS(h, x_poly, hasconst=True).fit()
    poly_h = np.poly1d(np.flipud(reg_h.params))

    replace_h = []
    for idx, row in dataframe_h0.iterrows():
        replace_h.append([idx, poly_h(np.log10(row['EMISSIE']))])

    if len(replace_h) > 0:
        replace_h = np.array(replace_h, dtype=float)
        if replace_h.ndim == 1:
            replace_h = replace_h.reshape(1, 2)

        if h_shape == 'log':
            replace_h[:, 1] = 10**replace_h[:, 1]
        elif h_shape == 'exp':
            replace_h[:, 1] = np.exp(replace_h[:, 1])

    del dataframe_h, dataframe_h0

    # ================================================================
    # 4. UITSTROOMOPENING_M2
    # ================================================================
    dataframe_a  = dataframe.loc[dataframe['UITSTROOMOPENING_M2'] > 0]
    dataframe_a0 = dataframe.loc[dataframe['UITSTROOMOPENING_M2'] == 0]

    replace_a = []

    # Only attempt regression if enough samples exist
    if len(dataframe_a) > 3:

        x = np.log10(dataframe_a.EMISSIE.values)

        if area_shape == 'log':
            a = np.log10(dataframe_a.UITSTROOMOPENING_M2.values)
        elif area_shape == 'exp':
            a = np.log(dataframe_a.UITSTROOMOPENING_M2.values)
        else:
            a = dataframe_a.UITSTROOMOPENING_M2.values

        x_poly = PolynomialFeatures(a_order).fit_transform(x.reshape(-1, 1))
        reg_a  = sm.OLS(a, x_poly).fit()
        poly_a = np.poly1d(np.flipud(reg_a.params))

        for idx, row in dataframe_a0.iterrows():
            replace_a.append([idx, poly_a(np.log10(row['EMISSIE']))])

    # Fallback if regression could not be used
    if len(replace_a) == 0 and len(dataframe_a0) > 0:

        # Engineering fallback using stack height
        for idx, row in dataframe_a0.iterrows():
            H = row['HOOGTE']
            if H > 0:
                # Diameter estimate: D = H / 30
                D = H / 30.0
                A = (np.pi / 4.0) * D * D
                replace_a.append([idx, A])

    # Convert safely to ndarray
    replace_a = np.array(replace_a, dtype=float)
    if replace_a.size > 0 and replace_a.ndim == 1:
        replace_a = replace_a.reshape(1, 2)

    # Undo transforms for regression-based fills
    if replace_a.size > 0 and len(dataframe_a) > 3:
        if area_shape == 'log':
            replace_a[:, 1] = 10**replace_a[:, 1]
        elif area_shape == 'exp':
            replace_a[:, 1] = np.exp(replace_a[:, 1])

    del dataframe_a, dataframe_a0

    # ================================================================
    return replace_v, replace_t, replace_h, replace_a



def write_df_to_csv(dataframe,datadir,year):
    dataframe.to_csv(datadir+f'point_source_gapfilled_{year}.csv')
    
    
    
def df2list(df, xmin, xmax, ymin, ymax, dx, dy, tracer_index):
    df_selection = df[(df['XCO_EMISSIEPUNT_HARM'] > xmin) &
                      (df['XCO_EMISSIEPUNT_HARM'] <= xmax) &
                      (df['YCO_EMISSIEPUNT_HARM'] > ymin) &
                      (df['YCO_EMISSIEPUNT_HARM'] <= ymax) &
                      (df['VOLUMESTROOM'] > 0) &
                      (df['TEMPERATUUR'] > 0) &
                      (df['HOOGTE'] > 0)].copy()  # Create a copy of the DataFrame

    df_selection['XCO_EMISSIEPUNT_HARM_INDEX'] = np.floor((df_selection['XCO_EMISSIEPUNT_HARM'].values - xmin) / dx) + 2
    df_selection['YCO_EMISSIEPUNT_HARM_INDEX'] = np.floor((df_selection['YCO_EMISSIEPUNT_HARM'].values - ymin) / dy) + 2
    df_selection['TRACER_INDEX'] = tracer_index

    data = df_selection[['XCO_EMISSIEPUNT_HARM', 'XCO_EMISSIEPUNT_HARM_INDEX',
                         'YCO_EMISSIEPUNT_HARM', 'YCO_EMISSIEPUNT_HARM_INDEX',
                         'EMISSIE', 'VOLUMESTROOM', 'TEMPERATUUR', 'HOOGTE', 'TRACER_INDEX', 'SNAP']].to_numpy()
    data = np.nan_to_num(data)

    return data



def data2netc_old(data, targetdir, nprocx = 999, nprocy = 999, tracer ='null', hour = 99, day = 99, month = 99, year = 9999):
    
    import netCDF4 as netc
    import datetime as datetime
    
    
    tfname  = f'{targetdir}pointsources.{year:04d}{month:02d}{day:02d}{hour:02d}.{tracer}.x{nprocx:03d}.y{nprocy:03d}.nc'
    
    do_proceed = True
    if do_proceed:
        # --- Create netCDF file
        tfobj = netc.Dataset(tfname, 'w')

        # --- Global attributes
        tfobj.title   = f"{tracer} Point sources for online plume rise calculations"
        tfobj.history = "Created: " + datetime.datetime.now().strftime("%d %b %Y")

        tfobj.description = f"Processed data from Emissieregistratie point source data for {tracer}, contains information on " \
                            "emission location in E-W direction (Rijksdriehoek coordinates, m) " \
                            "emission location in S-N direction (Rijksdriehoek coordinates, m) " \
                            f"emission strength Kilogram {tracer.lower()} per hour (kg hour-1)" \
                            "volume flux, (m3/s)" \
                            "exhaust temperature (K)" \
                            "and chimney height (m)"
                             
        tfobj.valid       = f"Valid in the hour starting at: {year:04d}-{month:02d}-{day:02d} {hour:02d}"
        tfobj.author      = 'A. Doyennel (VU)'
        tfobj.email       = 'a.doyennel@vu.nl'

        # -- Declaration of dimensions and variables
        dim_n     = tfobj.createDimension('n', len(data))
        dim_attrs = tfobj.createDimension('a', 6)

        var_attrs = tfobj.createVariable('attrs', 'f8', ('n','a' ))

        # -- Assigning values to variables
        var_attrs[:,:] = data

        # -- Variable attributes

        tfobj.close()    
        

        
def data2netc(data, targetdir, author, email, nprocx = 999, nprocy = 999, tracer ='null', minute = 99,hour = 99, day = 99, month = 99, year = 9999):
    
    import netCDF4 as netc
    import datetime as datetime
    
    tfname  = f'{targetdir}pointsources.{year:04d}{month:02d}{day:02d}{hour:02d}{minute:02d}.x{nprocx:03d}.y{nprocy:03d}.nc'
    
    do_proceed = True
    if do_proceed:
        # --- Create netCDF file
        tfobj = netc.Dataset(tfname, 'w')

        # --- Global attributes
        tfobj.title   = f"Point sources for online plume rise calculations"
        tfobj.history = "Created: " + datetime.datetime.now().strftime("%d %b %Y")

        tfobj.description = f"Processed data from Emissieregistratie point source data, contains information on " \
                            "emission location in E-W direction (HARMONIE coordinates, m) " \
                            "emission location in S-N direction (HARMONIE coordinates, m) " \
                            f"emission strength kilogram per hour (kg hour-1)" \
                            "volume flux, (m3/s)" \
                            "exhaust temperature (K)" \
                            "and chimney height (m)"
                             
        tfobj.valid       = f"Valid in the hour starting at: {year:04d}-{month:02d}-{day:02d} {hour:02d}"
        tfobj.author      = author
        tfobj.email       = email

        # -- Declaration of dimensions and variables
        dim_n     = tfobj.createDimension('n', len(data))

        var_x  = tfobj.createVariable('x',           'f4', ('n' ))
        var_xi = tfobj.createVariable('x_idx',       'f4', ('n' ))
        
        var_y  = tfobj.createVariable('y',           'f4', ('n' ))
        var_yi = tfobj.createVariable('y_idx',       'f4', ('n' ))
        
        var_z  = tfobj.createVariable('height',      'f4', ('n' ))
        var_e  = tfobj.createVariable('emission',    'f4', ('n' ))
        var_v  = tfobj.createVariable('volume',      'f4', ('n' ))
        var_t  = tfobj.createVariable('temperature', 'f4', ('n' ))
        var_s  = tfobj.createVariable('tracer_idx',  'f4', ('n' ))
        
        # -- Assigning values to variables
        var_x[:]  = data[:,0]
        var_xi[:] = data[:,1]
        
        var_y[:]  = data[:,2]
        var_yi[:] = data[:,3]
        
        var_e[:]  = data[:,4]
        var_v[:]  = data[:,5]
        var_t[:]  = data[:,6]+273.15
        var_z[:]  = data[:,7]
        var_s[:]  = data[:,8]

        # -- Variable attributes
        var_x.units  = 'm'
        var_xi.units = 'unitless'
        var_y.units  = 'm'
        var_yi.units = 'unitless'
        var_z.units  = 'm'
        var_e.units  = f'kg hour-1'
        var_v.units  = 'm3/s'
        var_t.units  = 'K'
        var_s.units  = 'unitless'
        
        tfobj.close() 


def df2list_new(df, xt_list, yt_list, tracer_name, domain_bounds=None):
    """
    Converts dataframe to point source array using grid-based indexing from les_grid.

    Parameters:
        df (DataFrame): Input emission dataframe
        les_grid (GridDales): Grid object with xt, yt arrays
        tracer_name (str): Name of the tracer (e.g. 'co2', 'ch4')
        domain_bounds (tuple): Optional spatial filter (xmin, xmax, ymin, ymax)

    Returns:
        np.ndarray: Array with columns [x, x_idx, y, y_idx, emission, volume, temp, height, stack_exit_area, tracer_name]
    """

    # Optional spatial bounding box
    if domain_bounds:
        xmin, xmax, ymin, ymax = domain_bounds
        df = df[(df['XCO_EMISSIEPUNT_HARM'] > xmin) &
                (df['XCO_EMISSIEPUNT_HARM'] <= xmax) &
                (df['YCO_EMISSIEPUNT_HARM'] > ymin) &
                (df['YCO_EMISSIEPUNT_HARM'] <= ymax)]

    # Apply basic physical constraints
    df_selection = df[
        (df['VOLUMESTROOM'] > 0) &
        (df['TEMPERATUUR'] > 0) &
        (df['HOOGTE'] > 0) &
        (df['UITSTROOMOPENING_M2'] > 0)
    ].copy()

    if df_selection.empty:
        return np.empty((0, 9))  # Return empty array with correct shape

    # Clip to grid bounds
    x_min, x_max = xt_list.min(), xt_list.max()
    y_min, y_max = yt_list.min(), yt_list.max()

    df_selection = df_selection[
        (df_selection['XCO_EMISSIEPUNT_HARM'] >= x_min) &
        (df_selection['XCO_EMISSIEPUNT_HARM'] <= x_max) &
        (df_selection['YCO_EMISSIEPUNT_HARM'] >= y_min) &
        (df_selection['YCO_EMISSIEPUNT_HARM'] <= y_max)
    ].copy()

    if df_selection.empty:
        return np.empty((0, 9))

    # Get emission coordinates
    x_emis = df_selection['XCO_EMISSIEPUNT_HARM'].values
    y_emis = df_selection['YCO_EMISSIEPUNT_HARM'].values

    # Find closest grid index (DALES expects indices starting at 2)
    x_idx = np.array([np.argmin(np.abs(xt_list - x)) for x in x_emis]) + 2
    y_idx = np.array([np.argmin(np.abs(yt_list - y)) for y in y_emis]) + 2

    df_selection['X_IDX'] = x_idx
    df_selection['Y_IDX'] = y_idx
    df_selection['TRACER_NAME'] = tracer_name

    # Build array
    data = df_selection[['XCO_EMISSIEPUNT_HARM', 'X_IDX',
                         'YCO_EMISSIEPUNT_HARM', 'Y_IDX',
                         'EMISSIE', 'VOLUMESTROOM', 'TEMPERATUUR', 'HOOGTE', 'UITSTROOMOPENING_M2', 'TRACER_NAME', 'SNAP']].to_numpy()

    return np.nan_to_num(data)



def data2netc_point(data, author, email, targetdir, tracer='null', minute=99, hour=99, day=99, month=99, year=9999):
    import netCDF4 as netc
    import datetime as datetime
    import os

    tracer = tracer.lower()  # Normalize tracer name

    tfname = os.path.join(
        targetdir,
        f'pointsources.{year:04d}{month:02d}{day:02d}{hour:02d}{minute:02d}.{tracer}.nc'
    )

    with netc.Dataset(tfname, 'w') as tfobj:
        # Global attributes
        tfobj.title = f"Point sources for online plume rise calculations"
        tfobj.history = "Created: " + datetime.datetime.now().strftime("%d %b %Y")
        tfobj.description = (
            f"Processed data from Emissieregistratie for tracer {tracer}. "
            "Contains emission location (m), emission strength (kg/h), volume flux (m³/s), "
            "temperature (K), chimney height (m), and stack exit area (m2)."
        )
        tfobj.valid = f"Valid for: {year:04d}-{month:02d}-{day:02d} {hour:02d}:00"
        tfobj.author = author
        tfobj.email = email

        # Dimensions
        n = len(data)
        tfobj.createDimension('n', n)

        # Variables
        tfobj.createVariable('x',           'f4', ('n'))[:] = data[:, 0]
        tfobj.createVariable('x_idx',       'i4', ('n'))[:] = data[:, 1]
        tfobj.createVariable('y',           'f4', ('n'))[:] = data[:, 2]
        tfobj.createVariable('y_idx',       'i4', ('n'))[:] = data[:, 3]
        tfobj.createVariable('height',      'f4', ('n'))[:] = data[:, 7]
        tfobj.createVariable('stack_exit_area',      'f4', ('n'))[:] = data[:, 8]
        tfobj.createVariable('emission',    'f4', ('n'))[:] = data[:, 4]
        tfobj.createVariable('volume',      'f4', ('n'))[:] = data[:, 5]
        tfobj.createVariable('temperature', 'f4', ('n'))[:] = data[:, 6] + 273.15

        tfobj.tracer_name = tracer

        # Units
        tfobj.variables['x'].units           = 'm'
        tfobj.variables['x_idx'].units       = 'grid index'
        tfobj.variables['y'].units           = 'm'
        tfobj.variables['y_idx'].units       = 'grid index'
        tfobj.variables['height'].units      = 'm'
        tfobj.variables['stack_exit_area'].units      = 'm-2'
        tfobj.variables['emission'].units    = f'kg hour-1'
        tfobj.variables['volume'].units      = 'm3/s'
        tfobj.variables['temperature'].units = 'K'



