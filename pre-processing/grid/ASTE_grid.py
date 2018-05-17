from lib_mitgcm_grid import *
import matplotlib.pylab as plt
import cartopy as cart

dirinputs='/Users/raphael/WORK/ASTE/input_files/'
diroutputs='/Users/raphael/WORK/ASTE/test_results/'

# define one grid object per facet
aste1 = mitgcm_grid(tile=1)
aste2 = mitgcm_grid(tile=2)
aste3 = mitgcm_grid(tile=3)
aste4 = mitgcm_grid(tile=4)
aste5 = mitgcm_grid(tile=5)

# read grid from binary files
aste1.load_grid_from_bin(gridfile=dirinputs + 'tileNTILE.mitgrid')
aste2.load_grid_from_bin(gridfile=dirinputs + 'tileNTILE.mitgrid')
aste3.load_grid_from_bin(gridfile=dirinputs + 'tileNTILE.mitgrid')
aste4.load_grid_from_bin(gridfile=dirinputs + 'tileNTILE.mitgrid')
aste5.load_grid_from_bin(gridfile=dirinputs + 'tileNTILE.mitgrid')

# test plot
plt.figure(figsize=[12,12])
m = plt.axes(projection=cart.crs.Orthographic(central_longitude=-45, central_latitude=60))
m.plot(aste1.grid.XC[:-1:20,:-1:20],aste1.grid.YC[:-1:20,:-1:20],'ro',alpha=0.5,transform=cart.crs.PlateCarree(),label='facet 1')
m.plot(aste2.grid.XC[:-1:20,:-1:20],aste2.grid.YC[:-1:20,:-1:20],'ko',alpha=0.5,transform=cart.crs.PlateCarree(),label='facet 2')
m.plot(aste3.grid.XC[:-1:20,:-1:20],aste3.grid.YC[:-1:20,:-1:20],'go',alpha=0.5,transform=cart.crs.PlateCarree(),label='facet 3')
m.plot(aste4.grid.XC[:-1:20,:-1:20],aste4.grid.YC[:-1:20,:-1:20],'bo',alpha=0.5,transform=cart.crs.PlateCarree(),label='facet 4')
m.plot(aste5.grid.XC[:-1:20,:-1:20],aste5.grid.YC[:-1:20,:-1:20],'yo',alpha=0.5,transform=cart.crs.PlateCarree(),label='facet 5')
m.coastlines()
m.add_feature(cart.feature.LAND, facecolor='0.75')
gl = m.gridlines(draw_labels=False)
plt.legend()
plt.show()

# read mask from previous run (we'd want a better solution)
aste1.infer_mask_from_output_file(diroutputs+'T.0000000003.data',precision='single',spval=0)
aste2.infer_mask_from_output_file(diroutputs+'T.0000000003.data',precision='single',spval=0)
aste3.infer_mask_from_output_file(diroutputs+'T.0000000003.data',precision='single',spval=0)
aste4.infer_mask_from_output_file(diroutputs+'T.0000000003.data',precision='single',spval=0)
aste5.infer_mask_from_output_file(diroutputs+'T.0000000003.data',precision='single',spval=0)

#write grids to netcdf
aste1.write_grid_to_nc('ASTE_FACET1.nc')
aste2.write_grid_to_nc('ASTE_FACET2.nc')
aste3.write_grid_to_nc('ASTE_FACET3.nc')
aste4.write_grid_to_nc('ASTE_FACET4.nc')
aste5.write_grid_to_nc('ASTE_FACET5.nc')
