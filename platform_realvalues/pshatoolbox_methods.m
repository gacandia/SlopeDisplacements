function ME = pshatoolbox_methods(group)

if nargin==1
    switch group
        case 1, load methods_GMM ME
        case 2, load methods_MR ME
        case 3, load methods_spatial ME
        case 4, load methods_spectral ME
        case 5, load methods_psda ME
        case 6, load methods_libs ME
    end
else

    % Generate matfiles with methods
    filepath  = fileparts(which('pshatoolbox_methods.m'));
    filename1 = fullfile(filepath,'methods_GMM');
    filename2 = fullfile(filepath,'methods_MR');
    filename3 = fullfile(filepath,'methods_spatial');
    filename4 = fullfile(filepath,'methods_spectral');
    filename5 = fullfile(filepath,'methods_psda');
    filename6 = fullfile(filepath,'methods_libs');


    %% GMM
    % IM = Intensity Measure Id.
    %    =  0  for PGA      (Peak Ground Acceleration)
    %    = -1  for PGV      (Peak Ground Velocity)
    %    = -2  for PGD      (Permanent Ground Displacement)
    %    = -3  for Duration (Duration)
    %    = -4  for CAV      (Cummulative Absolute Velocity)
    %    = -5  for AI       (Arias Intensity)
    %    = -6  for VGI      (Incremental Ground Velocity)
    %    = To    for Sa     (Pseudoacceleration    at period To , To>0)
    %    = To+i  for Sv     (Pseudo Velocity       at period To , To>0)
    %    = To+2i for Sd     (Spectral Displacement at period To , To>0)
    %    = To+3i for H/V    (H/V ratio             at period To , To>0)

    % Rmetric is a vector of the distance metrics used by the GMPE (1 or
    % if Rmetric(i)=1 means that the i-th metric will be computed, and 0
    % not be computed. The metris available are:
    %
    % Rmetric    = [Rrup Rhyp Rjb Repi Zseis Rx Ry0 zhyp ztor zbor zbot]
    % 1.-  rrup  = closest diatance from site to rupture plane
    % 2.-  rhyp  = distance from the site to the hypocenter
    % 3.-  rjb   = Joyner-Boore distance, i.e., closest distance from sit
    %              surface projection of rupture area
    % 4.-  repi  = epicentral distance
    % 5.-  rvolc = Rrup portion (km) within the volcanic arc
    % 6.-  rx    = Horizontal distance from top of rupture measured perpe
    %              to fault strike
    % 7.-  ry0   = Horizontal distance off the end of the rupture measure
    %              parallel to strike
    % 8.-  zhyp  = depth of hypocenter
    % 9.-  ztor  = Depth to top of coseismic rupture (km)
    % 10.- zbor  = Depth to the bottom of the rupture plane (km)
    % 11.- zbot  = The depth to the bottom of the seismogenic crust (km)


    ME(1:90,1) = createObj(6);
    % --------------------REGULAR MODELS (1-100)----------------------------------------------------------------------------
    i=1;  ME(i).id=i; ME(i).label = 'Youngs et al. 1997';           ME(i).handle = @Youngs1997;              ME(i).mech=[1 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.01 0.075 0.1 0.2 0.3 0.4 0.5 0.75 1 1.5 2 3];
    i=2;  ME(i).id=i; ME(i).label = 'Atkinson & Boore, 2003';       ME(i).handle = @AtkinsonBoore2003;       ME(i).mech=[1 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.01 0.04 0.1 0.2 0.4 1 2 1/0.33];
    i=3;  ME(i).id=i; ME(i).label = 'Zhao et al. 2006';             ME(i).handle = @Zhao2006;                ME(i).mech=[1 1 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.02 0.05 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.25 1.5 2 2.5 3 4 5];
    i=4;  ME(i).id=i; ME(i).label = 'McVerry et al. 2006';          ME(i).handle = @Mcverry2006;             ME(i).mech=[1 1 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.001 0.075 0.1 0.2 0.3 0.4 0.5 0.75 1 1.5 2 3];
    i=5;  ME(i).id=i; ME(i).label = 'Boroschek et al. 2012';        ME(i).handle = @ContrerasBoroschek2012;  ME(i).mech=[1 0 0]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.01 0.04 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 2];
    i=6;  ME(i).id=i; ME(i).label = 'Abrahamson et al. 2016';       ME(i).handle = @BCHydro2012;             ME(i).mech=[1 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];
    i=7;  ME(i).id=i; ME(i).label = 'Abrahamson et al. 2018';       ME(i).handle = @BCHydro2018;             ME(i).mech=[1 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 0 1 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];
    i=8;  ME(i).id=i; ME(i).label = 'Kuehn et al. 2020';            ME(i).handle = @Kuehn2020;               ME(i).mech=[1 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 0 1 0 0]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IM = [-1 0 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    i=9;  ME(i).id=i; ME(i).label = 'Parker et al. 2020';           ME(i).handle = @Parker2020;              ME(i).mech=[1 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IM = [-1 0.001 0.01 0.02 0.025 0.03 0.04 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 2.5 3 4 5 7.5 10];
    i=10; ME(i).id=i; ME(i).label = 'Abrahamson & Gulerce 2020';    ME(i).handle = @AG2020;                  ME(i).mech=[1 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 0 1 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.60 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];
    i=11; ME(i).id=i; ME(i).label = 'Idini et al. 2016';            ME(i).handle = @Idini2016;               ME(i).mech=[1 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.001 0.01 0.02 0.03 0.05 0.07 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    i=12; ME(i).id=i; ME(i).label = 'Montalva et al. 2017';         ME(i).handle = @MontalvaBastias2017;     ME(i).mech=[1 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];
    i=13; ME(i).id=i; ME(i).label = 'Montalva et al. 2017 (HQ) ';   ME(i).handle = @MontalvaBastias2017HQ;   ME(i).mech=[1 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];
    i=14; ME(i).id=i; ME(i).label = 'Montalva et al. 2018';         ME(i).handle = @Montalva2018;            ME(i).mech=[1 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];
    i=15; ME(i).id=i; ME(i).label = 'SIBER-RISK 2019';              ME(i).handle = @SiberRisk2019;           ME(i).mech=[1 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IM = [-1 0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];
    i=16; ME(i).id=i; ME(i).label = 'Garcia et al. 2005';           ME(i).handle = @Garcia2005;              ME(i).mech=[0 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IM = [-1 0 1/25 1/20 1/13.33 1/10 1/5 1/3.33 1/2.5 1/2 1/1.33 1/1 1/0.67 1/0.5 1/0.33 1/0.25 1/0.2];
    i=17; ME(i).id=i; ME(i).label = 'Jaimes et al. 2006';           ME(i).handle = @Jaimes2006;              ME(i).mech=[1 0 0]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.01 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 3.2 3.4 3.6 3.8 4 4.2 4.4 4.6 4.8 5 5.2 5.4 5.6 5.8 6];
    i=18; ME(i).id=i; ME(i).label = 'Jaimes et al. 2015';           ME(i).handle = @Jaimes2015;              ME(i).mech=[0 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IM = [-1 0.01 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 3.2 3.4 3.6 3.8 4 4.2 4.4 4.6 4.8 5];
    i=19; ME(i).id=i; ME(i).label = 'Jaimes et al. 2016';           ME(i).handle = @Jaimes2016;              ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IM = [-1 0.001 0.01 0.02 0.03 0.04 0.05 0.08 0.1 0.12 0.15 0.17 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    i=20; ME(i).id=i; ME(i).label = 'Garcia-Soto Jaimes 2017';      ME(i).handle = @GarciaJaimes2017;        ME(i).mech=[1 0 0]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IM = [-1 0.001 0.01 0.02 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3 4 5];
    i=21; ME(i).id=i; ME(i).label = 'Garcia-Soto Jaimes 2017 (V/H)';ME(i).handle = @GarciaJaimes2017VH;      ME(i).mech=[1 0 0]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[0 0 0 0 0 0 0 0 0 0 1]; ME(i).IM = 3j+[0.01 0.02 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3 4 5];
    i=22; ME(i).id=i; ME(i).label = 'Gulerce, Abrahamson 2011';     ME(i).handle = @GA2011;                  ME(i).mech=[1 1 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[0 0 0 0 0 0 0 0 0 0 1]; ME(i).IM = 3j+[-1 0.001 0.01 0.02 0.029 0.04 0.05 0.075 0.1 0.15 0.2 0.26 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    i=23; ME(i).id=i; ME(i).label = 'Stewart et al. 2016';          ME(i).handle = @SBSA2016;                ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[0 0 1 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IM = [-1,0.001,0.01,0.02,0.022,0.025,0.029,0.03,0.032,0.035,0.036,0.04,0.042,0.044,0.045,0.046,0.048,0.05,0.055,0.06,0.065,0.067,0.07,0.075,0.08,0.085,0.09,0.095,0.1,0.11,0.12,0.13,0.133,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.22,0.24,0.25,0.26,0.28,0.29,0.3,0.32,0.34,0.35,0.36,0.38,0.4,0.42,0.44,0.45,0.46,0.48,0.5,0.55,0.6,0.65,0.667,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.2,2.4,2.5,2.6,2.8,3,3.2,3.4,3.5,3.6,3.8,4,4.2,4.4,4.6,4.8,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10];
    i=24; ME(i).id=i; ME(i).label = 'Gulerce et al. 2017';          ME(i).handle = @GKAS2017;                ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 1 0 0 1 1 0 1 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 6 7.5 10];
    i=25; ME(i).id=i; ME(i).label = 'Bernal et al. 2014';           ME(i).handle = @Bernal2014;              ME(i).mech=[1 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.001 0.05 0.1 0.15 0.3 0.5 1 1.5 2 2.5 3 7 8];
    i=26; ME(i).id=i; ME(i).label = 'Sadigh et al. 1997';           ME(i).handle = @Sadigh1997;              ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.02 0.075 0.10 0.20 0.30 0.40 0.50 0.75 1 1.50 2 3 4];
    i=27; ME(i).id=i; ME(i).label = 'Idriss 2008';                  ME(i).handle = @I2008;                   ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    i=28; ME(i).id=i; ME(i).label = 'Chiou Youngs 2008';            ME(i).handle = @CY2008;                  ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 1 0 0 1 0 0 1 0 0]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IM = [-1 0.001 0.01 0.02 0.03 0.04 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    i=29; ME(i).id=i; ME(i).label = 'Boore Atkinson 2008';          ME(i).handle = @BA2008;                  ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[0 0 1 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IM = [-1 0.001 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    i=30; ME(i).id=i; ME(i).label = 'Campbell Bozorgnia 2008';      ME(i).handle = @CB2008;                  ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 1 0 0 0 0 0 1 0 0]; ME(i).IMtypes=[1 1 1 0 0 0 0 1 0 0 0]; ME(i).IM = [-2 -1 0.001 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    i=31; ME(i).id=i; ME(i).label = 'Abrahamson Silva 2008';        ME(i).handle = @AS2008;                  ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 1 0 0 1 0 0 1 0 0]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IM = [-1 0 0.01 0.02 0.03 0.04 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    i=32; ME(i).id=i; ME(i).label = 'Abrahamson Silva 1997 (Horz)'; ME(i).handle = @AS1997h;                 ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.01 0.02 0.03 0.04 0.05 0.06 0.075 0.09 0.1 0.12 0.15 0.17 0.2 0.24 0.3 0.36 0.4 0.46 0.5 0.6 0.75 0.85 1 1.5 2 3 4 5];
    i=33; ME(i).id=i; ME(i).label = 'Idriss 2014';                  ME(i).handle = @I2014;                   ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.001 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    i=34; ME(i).id=i; ME(i).label = 'CY 2014';                      ME(i).handle = @CY2014;                  ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 1 0 0 1 0 0 1 0 0]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IM = [-1 0.001 0.01 0.02 0.03 0.04 0.05 0.075 0.1 0.12 0.15 0.17 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    i=35; ME(i).id=i; ME(i).label = 'CB 2014';                      ME(i).handle = @CB2014;                  ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 1 0 0 1 0 1 1 0 1]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IM = [-1 0.001 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    i=36; ME(i).id=i; ME(i).label = 'BSSA 2014';                    ME(i).handle = @BSSA2014;                ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[0 0 1 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IM = [-1 0.001 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    i=37; ME(i).id=i; ME(i).label = 'ASK 2014';                     ME(i).handle = @ASK2014;                 ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 1 0 0 1 1 0 1 0 0]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IM = [-1 0.001 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 6 7.5 10];
    i=38; ME(i).id=i; ME(i).label = 'Akkar & Boomer 2007';          ME(i).handle = @AkkarBoomer2007;         ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[0 0 1 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 0 0 1 0]; ME(i).IM = [0 [0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2 2.05 2.1 2.15 2.2 2.25 2.3 2.35 2.4 2.45 2.5 2.55 2.6 2.65 2.7 2.75 2.8 2.85 2.9 2.95 3 3.05 3.1 3.15 3.2 3.25 3.3 3.35 3.4 3.45 3.5 3.55 3.6 3.65 3.7 3.75 3.8 3.85 3.9 3.95 4]+2*1j];
    i=39; ME(i).id=i; ME(i).label = 'Akkar & Boomer 2010';          ME(i).handle = @AkkarBoomer2010;         ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[0 0 1 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IM = [-1 0:0.05:3];
    i=40; ME(i).id=i; ME(i).label = 'Akkar et al. 2014';            ME(i).handle = @Akkar2014;               ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[0 1 1 1 0 0 0 0 0 0 0]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IM = [-1 0.001 0.01 0.02 0.03 0.04 0.05 0.075 0.1 0.15 0.2 0.3 0.4 0.5 0.75 1 1.5 2 3 4];
    i=41; ME(i).id=i; ME(i).label = 'Arroyo et al. 2010';           ME(i).handle = @Arroyo2010;              ME(i).mech=[1 0 0]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0 0.04 0.045 0.05 0.055 0.06 0.065 0.07 0.075 0.08 0.085 0.09 0.095 0.1 0.12 0.14 0.16 0.18 0.2 0.22 0.24 0.26 0.28 0.3 0.32 0.34 0.36 0.38 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.5 3 3.5 4 4.5 5];
    i=42; ME(i).id=i; ME(i).label = 'Bindi et al. 2011';            ME(i).handle = @Bindi2011;               ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[0 0 1 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IM = [-1 0 0.04 0.07 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8 0.9 1 1.25 1.5 1.75 2];
    i=43; ME(i).id=i; ME(i).label = 'Kanno et al. 2006';            ME(i).handle = @Kanno2006;               ME(i).mech=[1 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IM = [-1 0 0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.15 0.17 0.2 0.22 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.5 1.7 2 2.2 2.5 3 3.5 4 4.5 5];
    i=44; ME(i).id=i; ME(i).label = 'Cauzzi et al. 2015';           ME(i).handle = @Cauzzi2015;              ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 1 0 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[1 1 0 0 0 0 0 1 0 1 0]; ME(i).IM = [-1 round([0.01:0.01:0.09 0.1:0.05:10]*100)/100 round([0.01:0.01:0.09 0.1:0.05:10]*100)/100+2j];
    i=45; ME(i).id=i; ME(i).label = 'Du & Wang, 2012';              ME(i).handle = @DW12;                    ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[0 0 0 0 1 0 0 0 0 0 0]; ME(i).IM = -4;
    i=46; ME(i).id=i; ME(i).label = 'Foulser-Piggott, Goda 2015';   ME(i).handle = @FG15;                    ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[0 0 0 0 1 1 0 0 0 0 0]; ME(i).IM = [-5 -4];
    i=47; ME(i).id=i; ME(i).label = 'Travasarou et al. 2003';       ME(i).handle = @TBA03;                   ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[0 0 0 0 0 1 0 0 0 0 0]; ME(i).IM = -5;
    i=48; ME(i).id=i; ME(i).label = 'Bullock et al, 2017';          ME(i).handle = @BU17;                    ME(i).mech=[1 1 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[0 0 0 0 1 1 1 0 0 0 0]; ME(i).IM = [-6 -5 -4];
    i=49; ME(i).id=i; ME(i).label = 'Campbell,Bozorgnia 2010';      ME(i).handle = @CB10;                    ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 1 0 0 0 0 0 1 0 0]; ME(i).IMtypes=[0 0 0 0 1 0 0 0 0 0 0]; ME(i).IM = -4;
    i=50; ME(i).id=i; ME(i).label = 'Campbell,Bozorgnia 2011';      ME(i).handle = @CB11;                    ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 1 0 0 0 0 0 1 0 0]; ME(i).IMtypes=[0 0 0 0 1 0 0 0 0 0 0]; ME(i).IM = -4;
    i=51; ME(i).id=i; ME(i).label = 'Campbell,Bozorgnia 2019';      ME(i).handle = @CB19;                    ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 1 0 0 1 0 1 1 0 1]; ME(i).IMtypes=[0 0 0 0 1 1 0 0 0 0 0]; ME(i).IM = [-5 -4];
    i=52; ME(i).id=i; ME(i).label = 'Kramer & Mitchell, 2006';      ME(i).handle = @KM06;                    ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[0 0 0 0 1 0 0 0 0 0 0]; ME(i).IM = -4;
    i=53; ME(i).id=i; ME(i).label = 'PCE BCHydro (median)';         ME(i).handle = @medianPCEbchydro;        ME(i).mech=[1 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3];
    i=54; ME(i).id=i; ME(i).label = 'PCE NGA (median)';             ME(i).handle = @medianPCEnga;            ME(i).mech=[1 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.001 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 6 7 8 10];
    i=55; ME(i).id=i; ME(i).label = 'Jaimes et al. 2018';           ME(i).handle = @Jaimes2018;              ME(i).mech=[0 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.01 0.02 0.06 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 10];
    i=56; ME(i).id=i; ME(i).label = 'Arteta et al. 2018';           ME(i).handle = @Arteta2018;              ME(i).mech=[0 1 0]; ME(i).type=1;  ME(i).Rmetric=[0 1 0 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];
    i=57; ME(i).id=i; ME(i).label = 'Arteta 2023 (NoSAM-SC)';       ME(i).handle = @Arteta2023SC;            ME(i).mech=[0 0 1]; ME(i).type=1;  ME(i).Rmetric=[1 0 0 0 1 0 0 1 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 6 7.5 10];
    i=58; ME(i).id=i; ME(i).label = 'Arteta 2021 (NoSAM-SUB)';      ME(i).handle = @Arteta2021SUB;           ME(i).mech=[1 1 0]; ME(i).type=1;  ME(i).Rmetric=[1 1 0 0 0 0 0 1 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM = [0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 6 7.5 10];
    

    % --------------------PCE MODELS -----------------------------------------------------------------------
    i=70; ME(i).id=i; ME(i).label = 'PCE NGA';                      ME(i).handle = @PCE_nga;                 ME(i).mech=[0 0 1]; ME(i).type=2;  ME(i).Rmetric=[1 0 0 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM=[0.001 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 6 7 8 10];
    i=71; ME(i).id=i; ME(i).label = 'PCE BCHydro';                  ME(i).handle = @PCE_bchydro;             ME(i).mech=[1 1 0]; ME(i).type=2;  ME(i).Rmetric=[1 0 0 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[1 0 0 0 0 0 0 1 0 0 0]; ME(i).IM=[0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3];

    % --------------------SPECIAL MODELS  -----------------------------------------------------------------------
    i=75; ME(i).id=i; ME(i).label = 'Franky';                       ME(i).handle = @franky;                  ME(i).mech=[0 0 0]; ME(i).type=4;  ME(i).IM=[];
    i=76; ME(i).id=i; ME(i).label = 'Mean GMM';                     ME(i).handle = @meanGMM;                 ME(i).mech=[0 0 0]; ME(i).type=4;  ME(i).IM=[];

    % --------------------CONDITIONAL MODELS (200-)-----------------------------------------------------------------------
    % i=80; AI SUB
    % i=81; CAV shallow crustal
    % i=82; CAV SUB
    % i=83; CAVdp SUB
    i=80; ME(i).id=i; ME(i).label = 'Macedo, Abrahamson, Bray 2019';ME(i).handle = @MAB2019;                 ME(i).mech=[1 1 0]; ME(i).type=3;  ME(i).Rmetric=[1 0 0 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[0 0 0 0 0 1 0 0 0 0 0]; ME(i).IM=-5;
    i=81; ME(i).id=i; ME(i).label = 'Macedo, Abrahamson, Liu  2020';ME(i).handle = @MAL2020;                 ME(i).mech=[0 0 1]; ME(i).type=3;  ME(i).Rmetric=[1 0 1 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[0 0 0 0 1 0 0 0 0 0 0]; ME(i).IM=-4;
    i=82; ME(i).id=i; ME(i).label = 'Macedo, Liu, 2021';            ME(i).handle = @ML2021;                  ME(i).mech=[1 1 0]; ME(i).type=3;  ME(i).Rmetric=[1 0 0 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[0 0 0 0 1 0 0 0 0 0 0]; ME(i).IM=-4;
    i=83; ME(i).id=i; ME(i).label = 'Macedo et al., 2021';          ME(i).handle = @MCAVdp2021;              ME(i).mech=[1 1 0]; ME(i).type=3;  ME(i).Rmetric=[1 0 0 0 0 0 0 0 0 0 0]; ME(i).IMtypes=[0 0 0 0 1 0 0 0 0 0 0]; ME(i).IM=-4;

    % --------------------REGULAR MODELS (1-100)----------------------------------------------------------------------------
    i=1;  ME(i).ref='https://doi.org/10.1785/gssrl.68.1.58';
    i=2;  ME(i).ref='https://doi.org/10.1785/0120020156';
    i=3;  ME(i).ref='https://doi.org/10.1785/0120050122';
    i=4;  ME(i).ref='http://www.nzsee.org.nz/db/Bulletin/Archive/39(1)0001.pdf';
    i=5;  ME(i).ref='https://nisee.berkeley.edu/elibrary/files/documents/elib/www/documents/201204/PISELL/boroschek-maule-eq.pdf';
    i=6;  ME(i).ref='https://doi.org/10.1193/051712EQS188MR';
    i=7;  ME(i).ref='https://peer.berkeley.edu/sites/default/files/2018_02_abrahamson_9.10.18.pdf';
    i=8;  ME(i).ref='https://peer.berkeley.edu/sites/default/files/2020_02_final_3.17.2020.pdf';
    i=9;  ME(i).ref='https://peer.berkeley.edu/sites/default/files/2020_02_final_3.17.2020.pdf';
    i=10; ME(i).ref='';
    i=11; ME(i).ref='https://doi.org/10.1007/s10518-016-0050-1';
    i=12; ME(i).ref='https://doi.org/10.1785/0120160221';
    i=13; ME(i).ref='https://doi.org/10.1785/0120160221';
    i=14; ME(i).ref='';
    i=15; ME(i).ref='https://doi.org/10.1177/8755293019891723';
    i=16; ME(i).ref='https://doi.org/10.1785/0120050072';
    i=17; ME(i).ref='https://doi.org/10.1080/13632460609350622';
    i=18; ME(i).ref='https://doi.org/10.1080/13632469.2015.1025926';
    i=19; ME(i).ref='https://doi.org/10.1785/0120150283';
    i=20; ME(i).ref='https://doi.org/10.1785/0120160273';
    i=21; ME(i).ref='https://doi.org/10.1785/0120160273';
    i=22; ME(i).ref='https://doi.org/10.1193/1.3651317';
    i=23; ME(i).ref='https://doi.org/10.1193/072114eqs116m';
    i=24; ME(i).ref='http://dx.doi.org/10.1193/121814EQS213M';
    i=25; ME(i).ref='https://doi.org/10.13140/2.1.2693.6641';
    i=26; ME(i).ref='https://doi.org/10.1785/gssrl.68.1.180';
    i=27; ME(i).ref='https://doi.org/10.1193/1.2924362';
    i=28; ME(i).ref='https://doi.org/10.1193/1.2894832';
    i=29; ME(i).ref='https://doi.org/10.1193/1.2830434';
    i=30; ME(i).ref='https://doi.org/10.1193/1.2857546';
    i=31; ME(i).ref='https://doi.org/10.1193/1.2924360';
    i=32; ME(i).ref='https://doi.org/10.1785/gssrl.68.1.94';
    i=33; ME(i).ref='https://doi.org/10.1193/070613EQS195M';
    i=34; ME(i).ref='https://doi.org/10.1193/072813eqs219m';
    i=35; ME(i).ref='https://doi.org/10.1193/062913eqs175m';
    i=36; ME(i).ref='https://doi.org/10.1193/070113eqs184m';
    i=37; ME(i).ref='https://doi.org/10.1193/070913eqs198m';
    i=38; ME(i).ref='https://doi.org/10.1002/eqe.679';
    i=39; ME(i).ref='https://doi.org/10.1785/gssrl.81.2.195';
    i=40; ME(i).ref='https://doi.org/10.1007/s10518-013-9461-4';
    i=41; ME(i).ref='https://doi.org/10.1007/s10950-010-9200-0';
    i=42; ME(i).ref='https://doi.org/10.1007/s10518-011-9313-z';
    i=43; ME(i).ref='http://dx.doi.org/10.1785/0120050138';
    i=44; ME(i).ref='https://doi.org/10.1007/s10518-014-9685-y';
    i=45; ME(i).ref='https://doi.org/10.1002/eqe.2266';
    i=46; ME(i).ref='https://doi.org/10.1785/0120140316';
    i=47; ME(i).ref='https://doi.org/10.1002/eqe.270';
    i=48; ME(i).ref='https://doi.org/10.1785/0120160388';
    i=49; ME(i).ref='https://doi.org/10.1193/1.3457158';
    i=50; ME(i).ref='https://doi.org/10.1016/j.nucengdes.2011.04.020';
    i=51; ME(i).ref='https://doi.org/10.1193/090818EQS212M';
    i=52; ME(i).ref='https://doi.org/10.1193/1.2194970';
    i=53; ME(i).ref='https://doi.org/10.1016/j.enggeo.2020.105786';
    i=54; ME(i).ref='https://doi.org/10.1016/j.enggeo.2020.105786';
    i=55; ME(i).ref='';
    i=56; ME(i).ref='https://doi.org/10.1193/102116EQS176M';

    % --------------------PCE MODELS -----------------------------------------------------------------------
    i=70; ME(i).ref='https://doi.org/10.1016/j.enggeo.2020.105786';
    i=71; ME(i).ref='https://doi.org/10.1016/j.enggeo.2020.105786';

    % --------------------SPECIAL MODELS  -----------------------------------------------------------------------
    i=75; ME(i).ref='';
    i=76; ME(i).ref='';

    % --------------------CONDITIONAL MODELS (200-)-----------------------------------------------------------------------
    % i=80; AI SUB
    % i=81; CAV shallow crustal
    % i=82; CAV SUB
    % i=83; CAVdp SUB
    i=80; ME(i).ref='https://doi.org/10.1785/0120180297';
    i=81; ME(i).ref='https://doi.org/10.1785/0120190321';
    i=82; ME(i).ref='under review';
    i=83; ME(i).ref='under review';

    for i=1:numel(ME)
        if ~isnan(ME(i).id)
            ME(i).str=func2str(ME(i).handle);
        end
    end

    save(filename1,'ME')


    %% MR

    clear ME
    ME(1:4,1)=struct('id',[],'label',[],'func',[],'str',[],'ref',[]);
    i=1;ME(i).id=i; ME(i).label = 'Delta';                   ME(i).func = @delta;
    i=2;ME(i).id=i; ME(i).label = 'Truncated Exponential';   ME(i).func = @truncexp;
    i=3;ME(i).id=i; ME(i).label = 'Truncated Normal';        ME(i).func = @truncnorm;
    i=4;ME(i).id=i; ME(i).label = 'Characteristic';          ME(i).func = @yc1985;

    for i=1:numel(ME)
        if ~isempty(ME(i).id)
            ME(i).str=func2str(ME(i).func);
        end
    end

    save(filename2,'ME')


    %% spatial

    clear ME
    ME(1:4,1)=struct('id',[],'label',[],'func',[],'str',[],'param',[],'tooltip','','ref',[]);
    i=1;ME(i).id=i; ME(i).label = 'none';                              ME(i).func= @spa_none;    ME(i).param={''};        ME(i).ref='www.google.com';
    i=2;ME(i).id=i; ME(i).label = 'Jayaram N. and Baker J.W. (2009)';  ME(i).func= @spa_JB2009;  ME(i).param={'yes','no'};ME(i).ref='https://doi.org/10.1002/eqe.922'; ME(i).tooltip='Is Vs30 Clustered?';
    i=3;ME(i).id=i; ME(i).label = 'Loth, C., and Baker, J. W. (2013)'; ME(i).func= @spa_LB2013;  ME(i).param={''};        ME(i).ref='https://doi.org/10.1002/eqe.2212';
    i=4;ME(i).id=i; ME(i).label = 'Candia et al. (2019)';              ME(i).func= @spa_SR2019;  ME(i).param={''};        ME(i).ref='www.google.com';

    for i=1:numel(ME)
        if ~isempty(ME(i).id)
            ME(i).str=func2str(ME(i).func);
        end
    end

    save(filename3,'ME')


    %% spectral

    clear ME
    ME(1:13,1)=struct('id',[],'label',[],'func',[],'str',[],'param',[],'dependency',[],'tooltip','','ref',[]);                      % mechanism - magnitude - direction
    i=1; ME(i).id=i; ME(i).label = 'none';                     ME(i).func = @corr_none;                 ME(i).param={''}; ME(i).dependency = [0 0 0]; ME(i).ref='www.google.com';
    i=2; ME(i).id=i; ME(i).label = 'Baker & Cornell 2006';     ME(i).func = @corr_BakerCornell2006;     ME(i).param={'default','perpendicular'}; ME(i).dependency = [0 0 1];    ME(i).ref='https://doi.org/10.1785/0120050060';ME(i).tooltip='Ground motion component';
    i=3; ME(i).id=i; ME(i).label = 'Baker & Jayaram 2008';     ME(i).func = @corr_BakerJayaram2008;     ME(i).param={''}; ME(i).dependency = [0 0 0]; ME(i).ref='https://doi.org/10.1193/1.2857544';
    i=4; ME(i).id=i; ME(i).label = 'Jayaram et al. 2011';      ME(i).func = @corr_JayaramBaker2011;     ME(i).param={'crustal','interface','intraslab'}; ME(i).dependency = [1 0 0];    ME(i).ref='https://doi.org/10.12989/eas.2011.2.4.357';ME(i).tooltip='Earthquake mechanism';
    i=5; ME(i).id=i; ME(i).label = 'Cimellaro 2013';           ME(i).func = @corr_Cimellaro2013;        ME(i).param={'horizontal','vertical'}; ME(i).dependency = [0 0 1];    ME(i).ref='https://doi.org/10.1002/eqe.2248';ME(i).tooltip='Ground motion component';
    i=6; ME(i).id=i; ME(i).label = 'ASK2014 - NGA West2';      ME(i).func = @corr_Abrahamson2014;       ME(i).param={'within','between'}; ME(i).dependency = [0 0 0];    ME(i).ref='https://doi.org/10.1193/070913eqs198m';ME(i).tooltip='Residual';
    i=7; ME(i).id=i; ME(i).label = 'Abrahamanson et al. 2016'; ME(i).func = @corr_BCHhydro2016;         ME(i).param={''}; ME(i).dependency = [0 0 0];    ME(i).ref='https://doi.org/10.1193/051712EQS188MR';
    i=8; ME(i).id=i; ME(i).label = 'Baker & Bradley 2017';     ME(i).func = @corr_BakerBradley2017;     ME(i).param={''}; ME(i).dependency = [0 0 0];    ME(i).ref='https://doi.org/10.1193/060716EQS095M';
    i=9; ME(i).id=i; ME(i).label = 'Jaimes & Candia 2019';     ME(i).func = @corr_JaimesCandia2019;     ME(i).param={''}; ME(i).dependency = [0 0 0];    ME(i).ref='https://doi.org/10.1193/080918EQS200M';
    i=10;ME(i).id=i; ME(i).label = 'Candia et al. 2019';       ME(i).func = @corr_Candia2019;           ME(i).param={'interface','intraslab'}; ME(i).dependency = [1 0 0];    ME(i).ref='https://doi.org/10.1177/8755293019891723';ME(i).tooltip='Earthquake mechanism';
    i=11;ME(i).id=i; ME(i).label = 'Goda & Atkinson 2009';     ME(i).func = @corr_GodaAtkinson2009;     ME(i).param={''}; ME(i).dependency = [0 0 0];    ME(i).ref='https://doi.org/10.1785/0120090007';
    i=12;ME(i).id=i; ME(i).label = 'Akkar & Sandikkaya 2014';  ME(i).func = @corr_Akkar2014;            ME(i).param={''}; ME(i).dependency = [0 0 0];    ME(i).ref='https://doi.org/10.1007/s10518-013-9537-1';
    i=13;ME(i).id=i; ME(i).label = 'Campbell & Bozorgnia 2019';ME(i).func = @corr_CambpelBozorgnia2019; ME(i).param={''}; ME(i).dependency = [0 1 0];    ME(i).ref='https://doi.org/10.1193/090818EQS212M';

    for i=1:numel(ME)
        if ~isempty(ME(i).id)
            ME(i).str=func2str(ME(i).func);
        end
    end

    save(filename4,'ME')


    %% psda

    clear ME
    ME(1:0,1)=struct('id',[],'label',[],'isregular',[],'func',[],'str',[],'integrator',[],'primaryIM',[],'Safactor',[],'ref',[],'job',[]);

    % ------------------------- STANDARD MODELS ------------------------------------------------
    i=0;
    % subduction
    i=i+1; ME(i).id=i;ME(i).label = 'BMT 2017 Sa(M)';      ME(i).func = @BMT2017M;     ME(i).mechanism = 'subduction'; ME(i).integrator=1;  ME(i).primaryIM='SA(T=1.5Ts)';     ME(i).isregular=true;  ME(i).ref = 'https://doi.org/10.1061/(ASCE)GT.1943-5606.0001833';
    i=i+1; ME(i).id=i;ME(i).label = 'MC 2022F (V-M)';      ME(i).func = @MC2022MF;     ME(i).mechanism = 'subduction'; ME(i).integrator=7;  ME(i).primaryIM='PGV-SA(T=1.3Ts)'; ME(i).isregular=true;  ME(i).ref = '';
    i=i+1; ME(i).id=i;ME(i).label = 'MC 2022S (V-M)';      ME(i).func = @MC2022MS;     ME(i).mechanism = 'subduction'; ME(i).integrator=7;  ME(i).primaryIM='PGV-SA(T=1.3Ts)'; ME(i).isregular=true;  ME(i).ref = '';

    % shallow crustal (bray et al)
    i=i+1; ME(i).id=i;ME(i).label = 'BT 2007 Sa';          ME(i).func = @BT2007;       ME(i).mechanism = 'crustal'; ME(i).integrator=2;  ME(i).primaryIM='SA(T=1.5Ts)';     ME(i).isregular=true;  ME(i).ref = 'https://doi.org/10.1061/(ASCE)1090-0241(2007)133:4(381)';
    i=i+1; ME(i).id=i;ME(i).label = 'BT 2007 Sa(M)';       ME(i).func = @BT2007M;      ME(i).mechanism = 'crustal'; ME(i).integrator=1;  ME(i).primaryIM='SA(T=1.5Ts)';     ME(i).isregular=true;  ME(i).ref = 'https://doi.org/10.1061/(ASCE)1090-0241(2007)133:4(381)';
    i=i+1; ME(i).id=i;ME(i).label = 'BM 2019 NonNF (M)';   ME(i).func = @BM2019M;      ME(i).mechanism = 'crustal'; ME(i).integrator=1;  ME(i).primaryIM='SA(T=1.3Ts)';     ME(i).isregular=true;  ME(i).ref = 'https://www.ce.berkeley.edu/people/faculty/bray/research';

    % shallow crustal (other)
    i=i+1; ME(i).id=i;ME(i).label = 'Jibson  2007 (M)';    ME(i).func = @J07M;         ME(i).mechanism = 'crustal'; ME(i).integrator=1;  ME(i).primaryIM='PGA';           ME(i).isregular=true;  ME(i).ref = 'https://www.sciencedirect.com/science/article/pii/S0013795207000300?via%3Dihub';
    i=i+1; ME(i).id=i;ME(i).label = 'Jibson  2007 Ia';     ME(i).func = @J07Ia;        ME(i).mechanism = 'crustal'; ME(i).integrator=2;  ME(i).primaryIM='AI';            ME(i).isregular=true;  ME(i).ref = 'https://www.sciencedirect.com/science/article/pii/S0013795207000300?via%3Dihub';
    i=i+1; ME(i).id=i;ME(i).label = 'RA 2011';             ME(i).func = @RA2011;       ME(i).mechanism = 'crustal'; ME(i).integrator=4;  ME(i).primaryIM='PGV-PGA';       ME(i).isregular=true;  ME(i).ref = 'https://www.sciencedirect.com/science/article/pii/S0013795210002553';
    i=i+1; ME(i).id=i;ME(i).label = 'RS 2009 (Scalar-M)';  ME(i).func = @RS09M;        ME(i).mechanism = 'crustal'; ME(i).integrator=1;  ME(i).primaryIM='PGA';           ME(i).isregular=true;  ME(i).ref = 'http://www.nzsee.org.nz/db/Bulletin/Archive/42(1)0018.pdf';
    i=i+1; ME(i).id=i;ME(i).label = 'RS 2009 (Vector)';    ME(i).func = @RS09V;        ME(i).mechanism = 'crustal'; ME(i).integrator=3;  ME(i).primaryIM='PGV-PGA';       ME(i).isregular=true;  ME(i).ref = 'http://www.nzsee.org.nz/db/Bulletin/Archive/42(1)0018.pdf';
    i=i+1; ME(i).id=i;ME(i).label = 'AM 1988';             ME(i).func = @AM1988;       ME(i).mechanism = 'crustal'; ME(i).integrator=2;  ME(i).primaryIM='PGA';           ME(i).isregular=true;  ME(i).ref = 'https://doi.org/10.1002/eqe.4290160704';

    % ------------------------- PC MODELS ------------------------------------------------
    i=i+1; ME(i).id=i;ME(i).label = 'BMT 2017 (PCE-M)';       ME(i).func = @BMT2017_cdmM; ME(i).mechanism = 'subduction'; ME(i).integrator=6;  ME(i).primaryIM='SA(T=1.5Ts)';     ME(i).isregular=false; ME(i).ref = 'https://doi.org/10.1061/(ASCE)GT.1943-5606.0001833';
    i=i+1; ME(i).id=i;ME(i).label = 'BT 2007 (PCE-M)';        ME(i).func = @BT2007_cdmM;  ME(i).mechanism = 'crustal';    ME(i).integrator=6;  ME(i).primaryIM='SA(T=1.5Ts)';     ME(i).isregular=false; ME(i).ref = 'https://doi.org/10.1061/(ASCE)1090-0241(2007)133:4(381)';
    i=i+1; ME(i).id=i;ME(i).label = 'BT 2007 (PCE)';          ME(i).func = @BT2007_cdm;   ME(i).mechanism = 'crustal';    ME(i).integrator=5;  ME(i).primaryIM='SA(T=1.5Ts)';     ME(i).isregular=false; ME(i).ref = 'https://doi.org/10.1061/(ASCE)1090-0241(2007)133:4(381)';


    % ------------------------- MACHINE LEARNING BASED MODELS ------------------------------------------------
%     % single IM models
%     i=i+1; ME(i).id=i;ME(i).label = 'PLSR Interface(M)';     ME(i).func = @PLSR_F;       ME(i).mechanism = 'subductionF'; ME(i).integrator=1;  ME(i).job='plsr_1im_model_face.joblib'; ME(i).primaryIM='SA(T=1.3Ts)';     ME(i).isregular=true;  ME(i).ref = '';
%     i=i+1; ME(i).id=i;ME(i).label = 'PLSR Intraslab(M)';     ME(i).func = @PLSR_S;       ME(i).mechanism = 'subductionS'; ME(i).integrator=1;  ME(i).job='plsr_1im_model_slab.joblib'; ME(i).primaryIM='SA(T=1.3Ts)';     ME(i).isregular=true;  ME(i).ref = '';
% 
%     %two IM models
%     i=i+1; ME(i).id=i;ME(i).label = 'ML testing Sa(M)';      ME(i).func = @MLtestmdl;    ME(i).mechanism = 'subduction';  ME(i).integrator=7;  ME(i).job='ML_test_mdl.joblib';         ME(i).primaryIM='PGV-SA(T=1.3Ts)'; ME(i).isregular=true;  ME(i).ref = '';
% 
%     i=i+1; ME(i).id=i;ME(i).label = 'Ridge Interface(M)';    ME(i).func = @RidgeMF;      ME(i).mechanism = 'subductionF'; ME(i).integrator=7;  ME(i).job='poly_ridge_modelF.joblib';   ME(i).primaryIM='PGV-SA(T=1.3Ts)';     ME(i).isregular=true;  ME(i).ref = '';
%     i=i+1; ME(i).id=i;ME(i).label = 'Ridge Intraslab(M)';    ME(i).func = @RidgeMS;      ME(i).mechanism = 'subductionS'; ME(i).integrator=7;  ME(i).job='poly_ridge_modelS.joblib';   ME(i).primaryIM='PGV-SA(T=1.3Ts)';     ME(i).isregular=true;  ME(i).ref = '';
% 
%     % not working, weird trends
%     i=i+1; ME(i).id=i;ME(i).label = 'SVR Interface(M)';      ME(i).func = @SVRMF;        ME(i).mechanism = 'subductionF'; ME(i).integrator=7;  ME(i).job={'nystroem_transformF.joblib','nystroem_transform_xF.joblib','svr_nystroem_modelF.joblib'}; ME(i).primaryIM='PGV-SA(T=1.3Ts)'; ME(i).isregular=true;  ME(i).ref = '';
%     i=i+1; ME(i).id=i;ME(i).label = 'SVR Intraslab(M)';      ME(i).func = @SVRMS;        ME(i).mechanism = 'subductionS'; ME(i).integrator=7;  ME(i).job={'nystroem_transformS.joblib','nystroem_transform_xS.joblib','svr_nystroem_modelS.joblib'}; ME(i).primaryIM='PGV-SA(T=1.3Ts)'; ME(i).isregular=true;  ME(i).ref = '';
% 
%     % joblib files too heavy
%     i=i+1; ME(i).id=i;ME(i).label = 'RF Interface(M)';       ME(i).func = @RFMF;         ME(i).mechanism = 'subductionF'; ME(i).integrator=7;  ME(i).job='random_forest_modelF.joblib'; ME(i).primaryIM='PGV-SA(T=1.3Ts)';     ME(i).isregular=true;  ME(i).ref = '';
%     i=i+1; ME(i).id=i;ME(i).label = 'RF Intraslab(M)';       ME(i).func = @RFMS;         ME(i).mechanism = 'subductionS'; ME(i).integrator=7;  ME(i).job='random_forest_modelS.joblib'; ME(i).primaryIM='PGV-SA(T=1.3Ts)';     ME(i).isregular=true;  ME(i).ref = '';
% 
%     i=i+1; ME(i).id=i;ME(i).label = 'GBDT Interface(M)';     ME(i).func = @GBDTMF;       ME(i).mechanism = 'subductionF'; ME(i).integrator=7;  ME(i).job='boost_modelF.json'; ME(i).primaryIM='PGV-SA(T=1.3Ts)';     ME(i).isregular=true;  ME(i).ref = '';
%     i=i+1; ME(i).id=i;ME(i).label = 'GBDT Intraslab(M)';     ME(i).func = @GBDTMS;       ME(i).mechanism = 'subductionS'; ME(i).integrator=7;  ME(i).job='boost_modelS.json'; ME(i).primaryIM='PGV-SA(T=1.3Ts)';     ME(i).isregular=true;  ME(i).ref = '';

    % model is not yet developed
    %     i=i+1; ME(i).id=i;ME(i).label = 'ANN Interface(M)';     ME(i).func = @ANNMF;         ME(i).mechanism = 'subductionF'; ME(i).integrator=7;  ME(i).job={'ann_transform_xF.joblib','annF.pth'}; ME(i).primaryIM='PGV-SA(T=1.3Ts)'; ME(i).isregular=true;  ME(i).ref = '';
    %     i=i+1; ME(i).id=i;ME(i).label = 'ANN Intraslab(M)';     ME(i).func = @ANNMS;         ME(i).mechanism = 'subductionS'; ME(i).integrator=7;  ME(i).job={'ann_transform_xS.joblib','annS.pth'}; ME(i).primaryIM='PGV-SA(T=1.3Ts)'; ME(i).isregular=true;  ME(i).ref = '';

%     i=i+1; ME(i).id=i;ME(i).label = 'PCR Interface(M)';     ME(i).func = @PCRMF;         ME(i).mechanism = 'subductionF'; ME(i).integrator=7;  ME(i).job={'pcr_transform_xF.joblib','pca_transformF.joblib','pcr_modelF.joblib'}; ME(i).primaryIM='PGV-SA(T=1.3Ts)'; ME(i).isregular=true;  ME(i).ref = '';
%     i=i+1; ME(i).id=i;ME(i).label = 'PCR Intraslab(M)';     ME(i).func = @PCRMS;         ME(i).mechanism = 'subductionS'; ME(i).integrator=7;  ME(i).job={'pcr_transform_xS.joblib','pca_transformS.joblib','pcr_modelS.joblib'}; ME(i).primaryIM='PGV-SA(T=1.3Ts)'; ME(i).isregular=true;  ME(i).ref = '';

    % pc models
    i=i+1; ME(i).id=i;ME(i).label = 'PLSR Interface(PCE-M)'; ME(i).func = @PLSR_cdmMF;   ME(i).mechanism = 'subduction'; ME(i).integrator=6;  ME(i).job='plsr_1im_model_face.joblib'; ME(i).primaryIM='SA(T=1.3Ts)';     ME(i).isregular=false; ME(i).ref = '';

    for j=1:length(ME)
        ME(j).Safactor=str2IM(regexp(ME(j).primaryIM,'\-','split'));
        ME(j).Safactor=ME(j).Safactor(:)';
    end

    for i=1:numel(ME)
        if ~isempty(ME(i).id)
            ME(i).str=func2str(ME(i).func);
        end
    end

    save(filename5,'ME')


    %% libs

    clear ME
    ME(1:10,1)=struct('id',[],'label',[],'func',[],'str',[],'IM',[],'integrator',[],'ref',[]);

    % Bray & Macedo 2017 LIBS
    i=1;ME(i).id=i;ME(i).label = 'Bray & Macedo 2017 (Ds)';              ME(i).func = @BrayMacedo2017Ds;     ME(i).IM='PGA-SA(T=1)-CAV'; ME(i).integrator=3; ME(i).ref = 'https://doi.org/10.1016/j.soildyn.2017.08.026';
    i=2;ME(i).id=i;ME(i).label = 'Juang et al. 2013  (Dv)';              ME(i).func = @Juang2013Dv;          ME(i).IM='PGA';             ME(i).integrator=0; ME(i).ref = 'https://doi.org/10.1016/j.soildyn.2017.08.026';
    i=3;ME(i).id=i;ME(i).label = 'Bray & Macedo 2017 (LIBS)';            ME(i).func = @BrayMacedo2017;       ME(i).IM='PGA-SA(T=1)-CAV'; ME(i).integrator=3; ME(i).ref = '';
    i=4;ME(i).id=i;ME(i).label = 'Hutabarat 2020 (De)';                  ME(i).func = @Hutabarat2020De;      ME(i).IM='PGA';             ME(i).integrator=0; ME(i).ref = '';

    % Bullock 2018 LIBS
    i=5;ME(i).id=i;ME(i).label = 'Bullock 2018 (Ds)';                    ME(i).func = @Bullock2018Ds;        ME(i).IM='CAV';             ME(i).integrator=1; ME(i).ref = 'https://doi.org/10.1680/jgeot.17.P.174';
    i=6;ME(i).id=i;ME(i).label = 'Bullock 2018 (Dv)';                    ME(i).func = @Bullock2018Dv;        ME(i).IM='CAV';             ME(i).integrator=1; ME(i).ref = 'https://doi.org/10.1680/jgeot.17.P.174';
    i=7;ME(i).id=i;ME(i).label = 'Bullock 2018 (LIBS)';                  ME(i).func = @Bullock2018;          ME(i).IM='CAV';             ME(i).integrator=1; ME(i).ref = 'https://doi.org/10.1680/jgeot.17.P.174';

    % Bullock 2018 Tilting
    i=8;ME(i).id=i;ME(i).label = 'Bullock, 2018 (Tilt Empirical)';       ME(i).func = @Bullock2018e_tilt;         ME(i).IM='CAV';             ME(i).integrator=1; ME(i).ref = 'https://doi.org/10.1680/jgeot.17.P.174';
    i=9;ME(i).id=i;ME(i).label = 'Bullock, 2018 (Tilt Semi Empirical)';  ME(i).func = @Bullock2018s_tilt;         ME(i).IM='CAV-VGI';         ME(i).integrator=1; ME(i).ref = 'https://doi.org/10.1680/jgeot.17.P.174';

    % Null Model
    i=10;ME(i).id=i;ME(i).label = 'Null';                                 ME(i).func = @libs_null;                 ME(i).IM='PGA';             ME(i).integrator=1; ME(i).ref = '';

    for i=1:numel(ME)
        if ~isempty(ME(i).id)
            ME(i).str=func2str(ME(i).func);
        end
    end

    save(filename6,'ME')

    %% exit gracefully
    ME=[];
end


