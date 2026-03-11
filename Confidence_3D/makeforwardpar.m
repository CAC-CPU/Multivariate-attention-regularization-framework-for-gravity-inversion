function forwardpar = makeforwardpar()

% forwardpar.Ho: Geomagnetic induction in the unit of nT
% forwardpar.Inclination: Inclination of geomagnetic field
% forwardpar.Declination: Declination of geomagnetic field
% forwardpar.Azimuth: Set to be zero except the following cases:
  % 1.model magnetic vector/tensor data && 
  % 2.x is not the same as Easting (Since in the code we take postive x as
  % easting; postive y as northing and z doward)

% forwardpar.component: potential field components to be modeled:
  % gravity:              'gz'
  % gravity tensor:       'gxx', 'gzy', 'gxz', 'gyx', 'gyy', 'gyz', 'gzx', 'gzy', 'gzz'
  % magnetic total field: 'ht'
  % magnetic vector:      'hx', 'hy', 'hz'
  % magnetic tensor:      'hxx', 'hzy', 'hxz', 'hyx', 'hyy', 'hyz', 'hzx', 'hzy', 'hzz'


% eg.: forwardpar.component=str2mat('ht') indicates magneitc total field
% forwardpar.component=str2mat('hx','hy','ht'); indicates modeling magnetic
% vector component hx hy and magneitc total field ht;


% forwardpar.dispflag:
  % 1: display the forward modeling result.
  % 0: do not display the forward modeling result.


forwardpar.Ho=60000;
forwardpar.Inclination=90;
forwardpar.Declination=0;
forwardpar.Azimuth=0;
forwardpar.component=char('gz');% different components seperated by a comma
% forwardpar.component=char('gz');% different components seperated by a comma
forwardpar.noise=5;
forwardpar.obslevel = 0; %default level as 0 indicates observation on topography.
forwardpar.dispflag = 1;
clear x y z M;
save forwardpar.mat forwardpar;