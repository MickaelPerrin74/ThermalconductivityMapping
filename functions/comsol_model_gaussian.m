function T = comsol_model_gaussian(varargin)

pos_x = varargin{1};
pos_y= varargin{2};
Params = varargin{3};

%% import comsol
import com.comsol.model.*
import com.comsol.model.util.*

%% create model
ModelUtil.tags;
ModelUtil.remove('Model');
model = ModelUtil.create('Model');
model.modelPath(cd);

%% create parameters
ModelUtil.tags;
model.param.set('pos_x', sprintf('%1.2f [um]', pos_x), 'Laser beam position X');
model.param.set('pos_y', sprintf('%1.2f [um]', pos_y), 'Laser beam position Y');
model.param.set('R', sprintf('%1.2f [um]', Params.membrane_radius), 'Hole radius');
model.param.set('r0', sprintf('%1.2f [um]', Params.laser_spot_radius), 'Beam spot radius');
model.param.set('r_tot', sprintf('%1.2f [um]', Params.laser_spot_radius_tot), 'Total beam spot radius');
model.param.set('R_supp', sprintf('%1.2f [um]', Params.support_radius), 'for support area');
model.param.set('t_gr', '0.353 [nm]', 'Graphene thickness');
% model.param.set('Q_abs', sprintf('%1.2f [mW]', Absorption * Params.power), 'Absorbed laser power');
model.param.set('T0', sprintf('%1.2f [K]', Params.temperature), 'Ambient temperature');
model.param.set('h_conv', sprintf('%1.4e [W/(K*m^2)]', Params.convection), 'Convection coefficient graphene to air');
model.param.set('g', sprintf('%1.4e [W/(K*m^2)]', Params.g), 'Interface thermal conductance per unit area between graphene and Au/Si3N4 substrate');

%% create model
ModelUtil.tags;
model.component.create('mod1', false);

model.component('mod1').geom.create('geom1', 2);
model.component('mod1').label('Model 1');
model.component('mod1').defineLocalCoord(false);

%% create 2D kappa map
model.component('mod1').func.create('int1', 'Interpolation');
model.component('mod1').func('int1').set('source', 'file');
model.component('mod1').func('int1').set('importedname', 'kappa.txt');
model.component('mod1').func('int1').set('importedstruct', 'Spreadsheet');
model.component('mod1').func('int1').set('importeddim', '2D');
model.component('mod1').func('int1').set('funcs', {'kappa' '1'});
model.component('mod1').func('int1').set('defvars', true);
model.component('mod1').func('int1').set('interp', 'neighbor');
model.component('mod1').func('int1').set('argunit', 'um');
model.component('mod1').func('int1').set('filename', [pwd '/tmp/kappa.txt']);
model.component('mod1').func('int1').importData;
model.component('mod1').func('int1').set('nargs', '2');
model.component('mod1').func('int1').set('struct', 'spreadsheet');

%% create Absorption map
model.component('mod1').func.create('int2', 'Interpolation');
model.component('mod1').func('int2').set('source', 'file');
model.component('mod1').func('int2').set('importedname', 'absorption.txt');
model.component('mod1').func('int2').set('importedstruct', 'Spreadsheet');
model.component('mod1').func('int2').set('importeddim', '2D');
model.component('mod1').func('int2').set('funcs', {'Q_abs' '1'});
model.component('mod1').func('int2').set('defvars', true);
model.component('mod1').func('int2').set('interp', 'neighbor');
model.component('mod1').func('int2').set('argunit', 'um');
model.component('mod1').func('int2').set('filename', [pwd '/tmp/absorption.txt']);
model.component('mod1').func('int2').importData;
model.component('mod1').func('int2').set('nargs', '2');
model.component('mod1').func('int2').set('struct', 'spreadsheet');

%% create geometry
ModelUtil.tags;
geom1 = model.component('mod1').geom('geom1');

geom1.lengthUnit([native2unicode(hex2dec({'00' 'b5'}), 'unicode') 'm']);
geom1.repairTolType('relative');

geom1.create('Outer_domain', 'Circle');
geom1.feature('Outer_domain').set('pos', [0 0]);
geom1.feature('Outer_domain').set('r', 'R');

geom1.create('suspended_graphene', 'Circle');
geom1.feature('suspended_graphene').set('pos', [0 0]);
geom1.feature('suspended_graphene').set('r', 'R_supp');

geom1.create('laser_spot', 'Circle');
geom1.feature('laser_spot').set('pos', [pos_x pos_y]);
geom1.feature('laser_spot').set('r', 'r0');

geom1.create('laser_spot_total', 'Circle');
geom1.feature('laser_spot_total').set('pos', [pos_x pos_y]);
geom1.feature('laser_spot_total').set('r', 'r_tot');

geom1.run;

%% plot geom
% mphgeom(model,'geom1','facelabels','on')

%% identify domains
ModelUtil.tags;
N_domain = geom1.getNDomains;
upDown = geom1.getUpDown;

Area = zeros(N_domain,1);
for i = 1:N_domain
    model.component('mod1').geom('geom1').measureFinal.selection.geom('geom1', 2);
    model.component('mod1').geom('geom1').measureFinal.selection.set(i);
    Area(i) = geom1.measureFinal.getVolume;
end

[~ , idxs] = sort(Area);
index_graphene = idxs(end);
index_graphene_suspended = idxs(end-1);

%% create materials
ModelUtil.tags;
model.component('mod1').material.create('mat1', 'Common');
model.component('mod1').material('mat1').selection.set(1:N_domain);
model.component('mod1').material('mat1').label('Graphene');
model.component('mod1').material('mat1').propertyGroup('def').set('thermalconductivity', {'kappa(x,y)' '0' '0' '0' 'kappa(x,y)' '0' '0' '0' 'kappa(x,y)'});

%% create physics
ModelUtil.tags;
model.component('mod1').physics.create('ht', 'HeatTransfer', 'geom1');
model.component('mod1').physics('ht').create('temp1', 'TemperatureBoundary', 1);
outside_boundary_index = find(upDown(1,:)==0);              % check connections with void
model.component('mod1').physics('ht').prop('ShapeProperty').set('boundaryFlux_temperature', false);
model.component('mod1').physics('ht').prop('EquationForm').set('form', 'Stationary');
model.component('mod1').physics('ht').prop('ConsistentStabilization').set('heatCrosswindDiffusion', false);
model.component('mod1').physics('ht').prop('ConsistentStabilization').set('glim', '0.01[K]/ht.helem');
model.component('mod1').physics('ht').prop('ConsistentStabilization').set('StreamlineDiffusionOldForm', true);
%model.component('mod1').physics('ht').prop('RadiationProperty').set('fieldName', 'root.J');
model.component('mod1').physics('ht').prop('PhysicalModelProperty').set('outOfPlaneProperty', false);
model.component('mod1').physics('ht').prop('PhysicalModelProperty').set('dz', 1);

% set outer boundaries
model.component('mod1').physics('ht').feature('temp1').selection.set(outside_boundary_index);

% add heat sources
for i = 1:N_domain
    model.component('mod1').physics('ht').create(['hs' num2str(i)], 'HeatSource', 2);
end

% add convection and thermal loss to substrate
model.component('mod1').physics('ht').feature('hs1').set('Q0', '-(g * (T-T0) / t_gr) - (h_conv * (T-T0) / t_gr)');
model.component('mod1').physics('ht').feature('hs1').selection.set(index_graphene);

% add convection and thermal loss to suspended graphene (factor 2 due to both sides beig exposed)
model.component('mod1').physics('ht').feature('hs2').set('Q0', '-2*(h_conv * (T-T0) / t_gr)');
model.component('mod1').physics('ht').feature('hs2').selection.set(index_graphene_suspended);


% add laser for no intersecting domains
if N_domain == 4
    laser_area = pi * Params.laser_spot_radius * Params.laser_spot_radius;
    index_laser = find(min(min(abs(Area - laser_area))) == abs(Area - laser_area));

    index_laser_tot = find(not(ismember(1:N_domain,[index_graphene index_graphene_suspended index_laser])));

    model.component('mod1').physics('ht').feature('hs3').selection.set(index_laser);
    model.component('mod1').physics('ht').feature('hs4').selection.set(index_laser_tot);

    if sqrt(pos_x^2 + pos_y^2) < Params.membrane_radius
        model.component('mod1').physics('ht').feature('hs3').set('Q0', '(Q_abs / (2*pi*r0^2) / t_gr * exp(- ((x-pos_x)^2)/(2*r0^2) -((y-pos_y)^2)/(2*r0^2))) - (2 * h_conv * (T-T0) / t_gr)'); % when suspended, factor for both sides of membrane
        model.component('mod1').physics('ht').feature('hs4').set('Q0', '(Q_abs / (2*pi*r0^2) / t_gr * exp(- ((x-pos_x)^2)/(2*r0^2) -((y-pos_y)^2)/(2*r0^2))) - (2 * h_conv * (T-T0) / t_gr)'); % when suspended, factor for both sides of membrane
    else
        model.component('mod1').physics('ht').feature('hs3').set('Q0', '(Q_abs / (2*pi*r0^2) / t_gr * exp(- ((x-pos_x)^2)/(2*r0^2) -((y-pos_y)^2)/(2*r0^2))) - (g * (T-T0) / t_gr) - (h_conv * (T-T0) / t_gr)'); % on support
        model.component('mod1').physics('ht').feature('hs4').set('Q0', '(Q_abs / (2*pi*r0^2) / t_gr * exp(- ((x-pos_x)^2)/(2*r0^2) -((y-pos_y)^2)/(2*r0^2))) - (g * (T-T0) / t_gr) - (h_conv * (T-T0) / t_gr)'); % on support
    end
end


% add laser for only outer spot intersecting
if N_domain == 5

    laser_area = pi * Params.laser_spot_radius * Params.laser_spot_radius;
    index_laser = find(min(min(abs(Area - laser_area))) == abs(Area - laser_area));

    model.component('mod1').physics('ht').feature('hs3').selection.set(index_laser);

    index_laser_tot = find(not(ismember(1:N_domain,[index_graphene index_graphene_suspended index_laser])));

    if sqrt(pos_x^2 + pos_y^2) > Params.membrane_radius
        if sum(upDown(1,:)==index_graphene & upDown(2,:)==index_laser_tot(1)) > 0
            index_laser_tot_in = index_laser_tot(2);
            index_laser_tot_out = index_laser_tot(1);
        else
            index_laser_tot_in = index_laser_tot(1);
            index_laser_tot_out = index_laser_tot(2);
        end
        model.component('mod1').physics('ht').feature('hs4').selection.set(index_laser_tot_in);
        model.component('mod1').physics('ht').feature('hs5').selection.set(index_laser_tot_out);

        model.component('mod1').physics('ht').feature('hs3').set('Q0', '(Q_abs / (2*pi*r0^2) / t_gr * exp(- ((x-pos_x)^2)/(2*r0^2) -((y-pos_y)^2)/(2*r0^2))) - (g * (T-T0) / t_gr) - (h_conv * (T-T0) / t_gr)'); % when supported,
        model.component('mod1').physics('ht').feature('hs4').set('Q0', '(Q_abs / (2*pi*r0^2) / t_gr * exp(- ((x-pos_x)^2)/(2*r0^2) -((y-pos_y)^2)/(2*r0^2))) - (2 * h_conv * (T-T0) / t_gr)'); % when suspended, factor for both sides of membrane
        model.component('mod1').physics('ht').feature('hs5').set('Q0', '(Q_abs / (2*pi*r0^2) / t_gr * exp(- ((x-pos_x)^2)/(2*r0^2) -((y-pos_y)^2)/(2*r0^2))) - (g * (T-T0) / t_gr) - (h_conv * (T-T0) / t_gr)'); % when supported

    else

        if sum(upDown(1,:)==index_graphene & upDown(2,:)==index_laser_tot(1)) > 0
            index_laser_tot_in = index_laser_tot(2);
            index_laser_tot_out = index_laser_tot(1);
        else
            index_laser_tot_in = index_laser_tot(1);
            index_laser_tot_out = index_laser_tot(2);
        end
        model.component('mod1').physics('ht').feature('hs4').selection.set(index_laser_tot_in);
        model.component('mod1').physics('ht').feature('hs5').selection.set(index_laser_tot_out);

        model.component('mod1').physics('ht').feature('hs3').set('Q0', '(Q_abs / (2*pi*r0^2) / t_gr * exp(- ((x-pos_x)^2)/(2*r0^2) -((y-pos_y)^2)/(2*r0^2))) - (2 * h_conv * (T-T0) / t_gr)'); % when suspended, factor for both sides of membrane
        model.component('mod1').physics('ht').feature('hs4').set('Q0', '(Q_abs / (2*pi*r0^2) / t_gr * exp(- ((x-pos_x)^2)/(2*r0^2) -((y-pos_y)^2)/(2*r0^2))) - (2 * h_conv * (T-T0) / t_gr)'); % when suspended, factor for both sides of membrane
        model.component('mod1').physics('ht').feature('hs5').set('Q0', '(Q_abs / (2*pi*r0^2) / t_gr * exp(- ((x-pos_x)^2)/(2*r0^2) -((y-pos_y)^2)/(2*r0^2))) - (g * (T-T0) / t_gr) - (h_conv * (T-T0) / t_gr)'); % on support

    end
end

% add laser for both spots intersecting
if N_domain == 6

    Sum_area = Area+Area';
    Sum_area = Sum_area.*~eye(size(Sum_area));

    Outer_area = pi * Params.laser_spot_radius_tot * Params.laser_spot_radius_tot - pi * Params.laser_spot_radius * Params.laser_spot_radius;
    [col, row] = find(min(min(abs(Sum_area - Outer_area))) == abs(Sum_area - Outer_area));

    index_laser_tot = [col(1) row(1)];

    index_laser = find(not(ismember(1:N_domain,[index_graphene index_graphene_suspended index_laser_tot])));

    if sum(upDown(1,:)==index_graphene & upDown(2,:)==index_laser_tot(1)) > 0
        index_laser_tot_in = index_laser_tot(2);
        index_laser_tot_out = index_laser_tot(1);
    else
        index_laser_tot_in = index_laser_tot(1);
        index_laser_tot_out = index_laser_tot(2);
    end

    if sum(upDown(1,:)==index_laser_tot_out & upDown(2,:)==index_laser(1)) > 0
        index_laser_in = index_laser(2);
        index_laser_out = index_laser(1);
    else
        index_laser_in = index_laser(1);
        index_laser_out = index_laser(2);
    end

    model.component('mod1').physics('ht').feature('hs3').selection.set(index_laser_in);
    model.component('mod1').physics('ht').feature('hs4').selection.set(index_laser_out);
    model.component('mod1').physics('ht').feature('hs5').selection.set(index_laser_tot_in);
    model.component('mod1').physics('ht').feature('hs6').selection.set(index_laser_tot_out);


    model.component('mod1').physics('ht').feature('hs3').set('Q0', '(Q_abs / (2*pi*r0^2) / t_gr * exp(- ((x-pos_x)^2)/(2*r0^2) -((y-pos_y)^2)/(2*r0^2))) - (2 * h_conv * (T-T0) / t_gr)'); % when suspended, factor for both sides of membrane
    model.component('mod1').physics('ht').feature('hs4').set('Q0', '(Q_abs / (2*pi*r0^2) / t_gr * exp(- ((x-pos_x)^2)/(2*r0^2) -((y-pos_y)^2)/(2*r0^2))) - (g * (T-T0) / t_gr) - (h_conv * (T-T0) / t_gr)'); % on support
    model.component('mod1').physics('ht').feature('hs5').set('Q0', '(Q_abs / (2*pi*r0^2) / t_gr * exp(- ((x-pos_x)^2)/(2*r0^2) -((y-pos_y)^2)/(2*r0^2))) - (2 * h_conv * (T-T0) / t_gr)'); % when suspended, factor for both sides of membrane
    model.component('mod1').physics('ht').feature('hs6').set('Q0', '(Q_abs / (2*pi*r0^2) / t_gr * exp(- ((x-pos_x)^2)/(2*r0^2) -((y-pos_y)^2)/(2*r0^2))) - (g * (T-T0) / t_gr) - (h_conv * (T-T0) / t_gr)'); % on support

end

% set heat properties graphene
model.component('mod1').physics('ht').feature('solid1').label('Heat Transfer in Solids');
model.component('mod1').physics('ht').feature('solid1').set('rho_mat', 'userdef');
model.component('mod1').physics('ht').feature('solid1').set('Cp_mat', 'userdef');

% set initial conditions
model.component('mod1').physics('ht').feature('init1').set('Tinit', 'T0');
model.component('mod1').physics('ht').feature('init1').label('Initial Values');
model.component('mod1').physics('ht').feature('temp1').set('T0', 'T0');

% set thermal insulation (not applicable)
model.component('mod1').physics('ht').feature('ins1').label('Thermal Insulation');

%% create mesh
ModelUtil.tags;
model.component('mod1').mesh.create('mesh1');

model.component('mod1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('mod1').mesh('mesh1').feature('size').set('hauto', Params.grid_resolution);
model.component('mod1').mesh('mesh1').run;

%% create study
ModelUtil.tags;
model.study.create('std1');
model.study('std1').create('stat', 'Stationary');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').create('d1', 'Direct');
model.sol('sol1').feature('s1').feature.remove('fcDef');

model.sol('sol1').attach('std1');
model.sol('sol1').label('Solver 1');
model.sol('sol1').feature('v1').feature('mod1_T').label('mod1.T');
model.sol('sol1').feature('s1').set('stol', '0.0010');
model.sol('sol1').feature('s1').feature('dDef').set('ooc', false);
model.sol('sol1').feature('s1').feature('fc1').set('initstep', 0.01);
model.sol('sol1').feature('s1').feature('fc1').set('minstep', 1.0E-6);
model.sol('sol1').feature('s1').feature('fc1').set('maxiter', 50);
model.sol('sol1').feature('s1').feature('fc1').set('termonres', false);
model.sol('sol1').feature('s1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s1').feature('d1').set('ooc', false);
model.sol('sol1').feature('s1').feature('d1').set('pardmtsolve', false);
model.sol('sol1').runAll;

%% get average spot temperature
ModelUtil.tags;
tabel1 = model.result.table.create('tbl1', 'Table');
model.result.numerical.create('av1', 'AvSurface');
model.result.numerical('av1').set('data', 'dset1');
model.result.numerical('av1').selection.set(index_laser);
model.result.numerical('av1').set('probetag', 'none');

model.result.dataset('dset1').set('frametype', 'spatial');
model.result.numerical('av1').set('table', 'tbl1');
model.result.numerical('av1').set('dataseries', 'average');
model.result.numerical('av1').setResult;

T = tabel1.getReal;

return
