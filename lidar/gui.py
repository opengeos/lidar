
import os
import pkg_resources
import richdem as rd
from .filtering import MedianFilter, MeanFilter, GaussianFilter
from .filling import ExtractSinks
from .slicing import DelineateDepressions
from .mounts import DelineateMounts
# from lidar import *
import PySimpleGUI as sg

def main():
    # identify the sample data directory of the package
    package_name = 'lidar'
    data_dir = pkg_resources.resource_filename(package_name, 'data/')

    # use the sample dem. Change it to your own dem if needed
    in_dem = os.path.join(data_dir, 'dem.tif')
    # set output directory. By default, use the temp directory under user's home directory
    out_dir = os.path.join(os.path.expanduser("~"), "temp")

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    with sg.FlexForm('lidar package GUI') as form:
        form_rows = [
            [sg.Text('Level-set Method for Delineating Topographic Hierarchy', size=(50, 1), font=("Arial", 14), text_color='black')],
            [sg.Text('Select DEM:', font=("Arial", 14))],
            [sg.InputText(in_dem, size=(60, 1)), sg.FileBrowse()],
            [sg.Text('Delineation Mode:', font=("Arial", 14))],
            [sg.Radio('Depressions', "RADIO1", default=True), sg.Radio('Mounts', "RADIO1")],
            [sg.Text('DEM Filtering:', font=("Arial", 14))],
            [sg.Text('Select Filter:'), sg.InputCombo(['None', 'Mean Filter', 'Median Filter', 'Gaussian Filter']), sg.Text('Kernel Size: '), sg.InputText(default_text='3', size=(10, 1))],
            [sg.Text('Level-set Parameters:', font=("Arial", 14))],
            [sg.Text('Minimum size:'), sg.InputText(default_text='1000', size=(10, 1)), sg.Text('Minimum depth:'), sg.InputText(default_text='1.0', size=(10, 1))],
            [sg.Text('Slicing interval:'), sg.InputText(default_text='0.5', size=(10, 1)), sg.Text('Output shapefiles:'), sg.InputCombo(['Yes', 'No'], default_value='No')],
            [sg.Text('Display Results:', font=("Arial", 14))],
            [sg.InputCombo(['Yes', 'No'], default_value='No')],
            [sg.Text('Select Output Directory:', font=("Arial", 14))],
                    [sg.InputText(out_dir, size=(60, 1)), sg.FolderBrowse()],
                    [sg.Submit(), sg.Cancel()]]
        button, (in_dem, mode_dep, mode_mnt, filter_type, kernel_szie, min_size, min_depth, interval, bool_shp, display, out_dir) = form.LayoutAndRead(form_rows)
        # sg.Popup(button, source_filename)

        if button == 'Submit':

            kernel_szie = int(kernel_szie)
            min_size = int(min_size)
            min_depth = float(min_depth)
            interval = float(interval)
            if bool_shp == 'Yes':
                bool_shp = True
            else:
                bool_shp = False
            if display == 'Yes':
                display = True
            else:
                display = False
            if mode_mnt and in_dem == os.path.join(data_dir, 'dem.tif'):
                in_dem = os.path.join(data_dir, 'dsm.tif')

            out_dem_name = filter_type.split(" ")[0].lower() + '.tif'
            out_dem = os.path.join(out_dir, out_dem_name)

            sg.Popup('Please Wait!', 'The program is running! You will receive another message when it is done!')

            if filter_type == 'Mean Filter':                
                in_dem = MeanFilter(in_dem, kernel_size=kernel_szie, out_file=out_dem)
            elif filter_type == 'Median Filter':
                in_dem = MedianFilter(in_dem, kernel_size=kernel_szie, out_file=out_dem)
            elif filter_type == 'Gaussian Filter':
                in_dem = GaussianFilter(in_dem, sigma=kernel_szie, out_file=out_dem)

            if mode_dep:
                sink_path = ExtractSinks(in_dem, min_size, out_dir)
                dep_id_path, dep_level_path = DelineateDepressions(sink_path, min_size, min_depth, interval, out_dir, bool_shp)
            else:
                sink_path = os.path.join(out_dir, 'sink.tif')
                dep_id_path, dep_level_path = DelineateMounts(in_dem, min_size, min_depth, interval, out_dir, bool_shp)

            if display:
                # loading data and results
                dem = rd.LoadGDAL(in_dem)
                sink = rd.LoadGDAL(sink_path)
                dep_id = rd.LoadGDAL(dep_id_path)
                dep_level = rd.LoadGDAL(dep_level_path)

                # plotting results
                dem_fig = rd.rdShow(dem, ignore_colours=[0], axes=False, cmap='jet', figsize=(6, 5.5))
                sink_fig = rd.rdShow(sink, ignore_colours=[0], axes=False, cmap='jet', figsize=(6, 5.5))
                dep_id_fig = rd.rdShow(dep_id, ignore_colours=[0], axes=False, cmap='jet', figsize=(6, 5.5))
                dep_level_path = rd.rdShow(dep_level, ignore_colours=[0], axes=False, cmap='jet', figsize=(6, 5.5))

                del dem, sink, dep_id, dep_level, dem_fig, sink_fig, dep_id_fig, dep_id_path

            sg.Popup('Success!','The results are saved in: {}'.format(out_dir))

def GUI():
    main()

def gui():
    main()

if __name__ == '__main__':
    main()