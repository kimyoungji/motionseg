function energy = GCO_GetEnergy(Handle,varargin)
% GCO_GetLabeling     Retrieve the current energy
%     GCO_GetLabeling(Handle,i,j) returns energy of site i when having
%     label l

GCO_LoadLib();
d_num = int32(varargin{1});
g_num = int32(varargin{2});
energy = gco_matlab('gco_getenergy',Handle,d_num,g_num);
end