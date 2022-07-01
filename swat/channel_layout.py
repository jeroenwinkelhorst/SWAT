import ipywidgets as widgets
from ipywidgets.widgets import HBox, VBox, Layout
from ipyfilechooser import FileChooser
from swat.channel import Trapezoidal
from swat.channel import YZ
import plotly.graph_objs as go
import pandas as pd
import numpy as np
import csv

def get_pos_trace(fig, legend_name):
    """Return position of trace in fig.data based on legend name"""
    for i in range(len(fig.data)):
        if fig.data[i].name == legend_name:
            return i


class ChannelNotebookLayout:
    def __init__(self):
        self.output = widgets.Output()
        self.geometry = widgets.Output()
        # widgets dict
        self.wg = {}
        self.wg_style = {}

        self.fig_dwp = None
        self.figure_geometry_legend = {}

        self._set_general_vars()
        self._widget_defs()
        self._create_fig_geometry()

    def __repr__(self):
        return "Channel notebook layout"

    def _set_general_vars(self):
        self.wg_style['fixed_width'] = {'description_width': '250px'}
        self.wg_style['layout'] = {'width': '350px'}

    def run(self):
        return VBox([
            widgets.HTML(value="<p style='font-family:verdana; font-size:24px; padding:20px'>Initial parameters</p>"),
            self.wg['cs_type'],
            self.output])

    def _add_general_fields(self):
        with self.output:
            display(
                HBox([VBox([
                    self.wg['h_or_q'],
                    self.wg['water_level'],
                    self.wg['discharge'],
                    self.wg['b_b'],
                    self.wg['k_manning'],
                    self.wg['bottom_slope'],
                    self.wg['surface_level'],
                    self.geometry
                ]),self.figdwp])
            )

    def _add_trapezoidal_fields(self):
        with self.geometry:
            display(
                VBox([
                    widgets.HTML(
                        value="<p style='font-family:verdana; font-size:24px; padding:20px'>Trapezoidal geometry</p>"),
                    self.wg['slope_left'],
                    self.wg['slope_right'],
                    self.wg['bed_width'],
                    self.wg['bed_level']
                ])
            )

    def _add_yz_fields(self):
        with self.geometry:
            display(
                VBox([
                    widgets.HTML(
                        value="<p style='font-family:verdana; font-size:24px; padding:20px'>YZ cross section geometry</p>"),
                    self.wg['yz_message'],    
                    self.wg['yz_file'],
                    self.wg['yz_but'],
                    self.out_yz
                ])
            )
                    

    def _add_result_fields(self):
        with self.output:
            display(
                VBox([widgets.HTML(value="<p style='font-family:verdana; font-size:24px; padding:20px'>Output</p>"),
                HBox([
                    VBox([
                        
                        self.wg['discharge_out'],
                        self.wg['water_level_out'],
                        self.wg['velocity'],
                        self.wg['wet_area'],
                        self.wg['wet_perimeter'],
                        self.wg['hydraulic_radius'],
                        self.wg['froude'],
                        self.wg['freeboard_left'],
                        self.wg['freeboard_right']
                    ]),
                    VBox([self.wg['warning'],self.wg['froude_message']])
                ])
                ])
            )

    def _observe_cs_type(self, b):
        self.wg['cs_type'].disabled = True
        self._add_general_fields()
        if self.wg['cs_type'].value == "Trapezoidal":
            self._add_trapezoidal_fields()
        if self.wg['cs_type'].value == "YZ":
            self._add_yz_fields()    
        self._add_result_fields()

    def _observe_cs_h_or_q(self, b):
        if self.wg['h_or_q'].value==True:
            self.wg['water_level'].disabled=True
            self.wg['discharge'].disabled=False
        else:
            self.wg['water_level'].disabled=False
            self.wg['discharge'].disabled=True      

    def _observe_cs_input(self, b):
        self._update()

    def _widget_defs(self):
        style = self.wg_style['fixed_width']
        layout =  self.wg_style['layout']

        self.wg['cs_type'] = widgets.Dropdown(
            options=['', 'Trapezoidal','YZ'],
            value='',
            description='Select cross section type:',
            disabled=False,
            style=style,
            layout = layout
        )

        self.wg['cs_type'].observe(self._observe_cs_type, names='value')

        self.wg['h_or_q'] = widgets.Checkbox(
            value=False,
            description='calculate water level',
            disabled=False,
            indent=False
            )
        self.wg['h_or_q'].observe(self._observe_cs_h_or_q, names='value') 
        self.wg['b_b'] = widgets.Checkbox(
            value=False,
            description='Use Bos en Bijkerk',
            disabled=False,
            indent=False
            )
        self.wg['b_b'].observe(self._observe_cs_input, names='value')    
        
        # General input
        self.wg['k_manning'] = widgets.BoundedFloatText(
            value=30,
            min=1,
            max=100,
            step=1,
            description='K strikler [m(1/3)/s]:',
            disabled=False,
            style=style,
            layout = layout
        )
        self.wg['k_manning'].observe(self._observe_cs_input, names='value')

        self.wg['bottom_slope'] = widgets.BoundedFloatText(
            value=0.001,
            min=0.0001,
            max=0.9999,
            step=0.0001,
            description='Channel slope [m/m]:',
            disabled=False,
            style=style,
            layout = layout
        )
        self.wg['bottom_slope'].observe(self._observe_cs_input, names='value')

        self.wg['water_level'] = widgets.BoundedFloatText(
            value=0.5,
            min=-10,
            max=2000,
            step=0.001,
            description='Water level [m+ref]:',
            disabled=False,
            style=style,
            layout = layout
        )
        self.wg['water_level'].observe(self._observe_cs_input, names='value')
        
        self.wg['discharge'] = widgets.BoundedFloatText(
            value=0.326,
            min=0,
            max=5000,
            step=0.001,
            description='Discharge [m3/s]:',
            disabled=True,
            style=style,
            layout = layout
        )
        
        self.wg['discharge'].observe(self._observe_cs_input, names='value')

        self.wg['surface_level'] = widgets.BoundedFloatText(
            value=1,
            min=-10,
            max=2000,
            step=0.01,
            description='Surface level [m+ref]:',
            disabled=False,
            style=style,
            layout = layout
        )
        self.wg['surface_level'].observe(self._observe_cs_input, names='value')

        # ouput

        self.wg['discharge_out'] = widgets.BoundedFloatText(
            value=0.326,
            max=2000000,
            description='Discharge [m3/s]:',
            disabled=True,
            style=style,
            layout = layout
        )

        self.wg['water_level_out'] = widgets.BoundedFloatText(
            value=0.5,
            max=2000,
            description='Water level [m+ref]:',
            disabled=True,
            style=style,
            layout = layout
        )

        self.wg['velocity'] = widgets.BoundedFloatText(
            value=0.326/0.75,
            min=-100,
            max=100,
            description='Velocity [m/s]:',
            disabled=True,
            style=style,
            layout = layout
        )

        self.wg['wet_area'] = widgets.BoundedFloatText(
            value=0.75,
            max=2000000,            
            description='Wet area [m2]:',
            disabled=True,
            style=style,
            layout = layout
        )
        #self.wg['wet_area'].observe(self._observe_cs_input, names='value')

        self.wg['wet_perimeter'] = widgets.BoundedFloatText(
            value=2.41,
            max=5000,
            
            description='Wet perimeter [m]:',
            disabled=True,
            style=style,
            layout = layout
        )
        
        self.wg['hydraulic_radius'] = widgets.BoundedFloatText(
            value=0.311,
            max=5000,
            
            description='Hydraulic radius [m]:',
            disabled=True,
            style=style,
            layout = layout
        )
        #self.wg['hydraulic_radius'].observe(self._observe_cs_input, names='value')

        self.wg['froude'] = widgets.BoundedFloatText(
            value=0.1965,
            max=5000,
            min=0,
            description='Froude number (1 = critical) [-]:',
            disabled=True,
            style=style,
            layout = layout
        )

        self.wg['freeboard_left'] = widgets.BoundedFloatText(
            value=0.5,
            max=500,
            min=-500,
            description='Freeboard (>0) or flood (<0) level left bank[m]:',
            disabled=True,
            style=style,
            layout = layout
        )

        self.wg['freeboard_right'] = widgets.BoundedFloatText(
            value=0.5,
            max=500,
            min=-500,
            description='Freeboard (>0) or flood (<0) level right bank[m]:',
            disabled=True,
            style=style,
            layout = layout
        )

        self.wg['warning']=widgets.HTML(value="",layout=Layout(display="flex", justify_content="flex-start", width="500px"))

        self.wg['froude_message']=widgets.HTML(value="",layout=Layout(display="flex", justify_content="flex-start", width="500px"))

        # Trapezoidal
        self.wg['slope_left'] = widgets.BoundedFloatText(
            value=1,
            min=0,
            max=100,
            step=0.01,
            description='slope left 1:n [n]:',
            disabled=False,
            style=style,
            layout = layout
        )
        self.wg['slope_left'].observe(self._observe_cs_input, names='value')

        self.wg['slope_right'] = widgets.BoundedFloatText(
            value=1,
            min=0,
            max=200,
            step=0.01,
            description='slope right 1:n [n]:',
            disabled=False,
            style=style,
            layout = layout
        )
        self.wg['slope_right'].observe(self._observe_cs_input, names='value')

        self.wg['bed_width'] = widgets.BoundedFloatText(
            value=1,
            min=0,
            max=200,
            step=0.01,
            description='bed width [m]:',
            disabled=False,
            style=style,
            layout = layout
        )
        self.wg['bed_width'].observe(self._observe_cs_input, names='value')

        self.wg['bed_level'] = widgets.BoundedFloatText(
            value=0,
            min=-10,
            max=2000,
            step=0.01,
            description='bed level [m]:',
            disabled=False,
            style=style,
            layout = layout
        )
        self.wg['bed_level'].observe(self._observe_cs_input, names='value')

        # YZ
        self.wg['yz_file']= FileChooser()
        self.wg['yz_but']=  widgets.Button(
            value=True,
            description='GO',
            button_style='warning',
            tooltip='Description',
            icon='fire')
        self.wg['yz_but'].on_click(self._observe_cs_input)    
        #self.wg['yz_file'].observe(self._observe_cs_input, names='value')

        self.out_yz = widgets.Output()

        self.wg['yz_message']=widgets.HTML(value="<p style='font-family:verdana; font-size:14px'>upload a .CSV files with 2 columns (Y and Z values, see YZ_example1.csv in example_files).</p>")

        
        

    def _create_fig_geometry(self):
        # Figure widget
        fig = go.FigureWidget()
        fig.layout.template = "plotly_white"
        fig.layout.xaxis.dtick = 1
        fig.update_xaxes()
        fig.update_xaxes(title_text='[m]')
        fig.update_yaxes(title_text='[m]')
        fig.update_yaxes(
            scaleanchor="x",
            scaleratio=1,
        )



        # channel
        legend_channel = "channel"
        fig.add_scatter(x=[0, 1, 2, 3], y=[1, 0, 0, 1],
                        line=dict(color="green"),
                        name=legend_channel)

        legend_wet_profile = "wet profile"
        fig.add_scatter(x=[0.5, 1, 2, 2.5], y=[0.5, 0, 0, 0.5],
                        line=dict(color="red"),
                        name=legend_wet_profile)
                # Water level
        legend_water_level = "water level"
        fig.add_scatter(x=[0, 3], y=[0.5, 0.5],
                        line=dict(color="blue"),
                        name=legend_water_level)                
        self.figure_geometry_legend['water level'] = legend_water_level
        self.figure_geometry_legend['channel'] = legend_channel
        self.figure_geometry_legend['wet_profile'] = legend_wet_profile
        self.figdwp = fig

    def _update(self):
        if self.wg['cs_type'].value == 'Trapezoidal':
            cs = Trapezoidal(
                k_manning = self.wg['k_manning'].value,
                b_b = self.wg['b_b'].value,
                bottom_slope = self.wg['bottom_slope'].value,
                slope_left = self.wg['slope_left'].value,
                slope_right = self.wg['slope_right'].value,
                bed_width = self.wg['bed_width'].value,
                bed_level = self.wg['bed_level'].value,
                surface_level = self.wg['surface_level'].value)
            if self.wg['h_or_q'].value == True:
                self.wg['discharge_out'].value = self.wg['discharge'].value
                self.wg['water_level_out'].value = cs.calc_water_level(self.wg['discharge'].value)
                self.wg['water_level'].value = cs.calc_water_level(self.wg['discharge'].value)
            else:    
                self.wg['discharge'].value = cs.calc_discharge(self.wg['water_level'].value)
                self.wg['discharge_out'].value = cs.calc_discharge(self.wg['water_level'].value)
                self.wg['water_level_out'].value = self.wg['water_level'].value
            self.wg['velocity'].value = cs.velocity
            self.wg['wet_area'].value = cs.wet_area
            self.wg['wet_perimeter'].value = cs.wet_perimeter
            self.wg['hydraulic_radius'].value = cs.hydraulic_radius
            self.wg['froude'].value = cs.froude_number
            self.wg['freeboard_left'].value = cs.freeboard
            self.wg['freeboard_right'].value = cs.freeboard
            if self.wg['water_level'].value<=self.wg['surface_level'].value:
                self.wg['warning'].value = ""
            else:

                self.wg['warning'].value = "<p style='font-family:verdana; font-size:16px'>Water level is larger than Surface level. Sides of part above surface level are taken as vertical (for area calculation). Perimeter and Hydraulic radius are based on wet profile (red line in graph). Make the surface level higher than the water level to correct</p>"
                
            self._update_fig_geometry(cs)
        elif self.wg['cs_type'].value == 'YZ':
            if self.wg['yz_file'].selected is not None:
                with open(self.wg['yz_file'].selected, 'rb') as csvfile:
                    
                    has_header = csv.Sniffer().has_header(csvfile.read(1024).decode('utf-8'))
                    if has_header:    
                        yz_table = pd.read_csv(self.wg['yz_file'].selected)
                    else:
                        yz_table = pd.read_csv(self.wg['yz_file'].selected,header=None)
                        yz_table.columns = {'Y', 'Z'}
                    self.out_yz.clear_output()
                    with self.out_yz:
                        display(widgets.HTML("Data from file:"))

                        display(yz_table)    
                    cs = YZ(
                        k_manning = self.wg['k_manning'].value,
                        b_b = self.wg['b_b'].value,
                        bottom_slope = self.wg['bottom_slope'].value,
                        yz_values = yz_table.values)
                    if cs.yz_ok():
                        self.wg['warning'].value = ""  
                        self.wg['yz_message'].value = "<p style='font-family:verdana; font-size:14px'>upload a .CSV files with 2 columns (Y and Z values, see YZ_example1.csv in example_files).</p>"
                        if self.wg['h_or_q'].value == True:
                            self.wg['discharge_out'].value = self.wg['discharge'].value
                            self.wg['water_level_out'].value = cs.calc_water_level(self.wg['discharge'].value)

                        else:    
                            self.wg['discharge_out'].value = cs.calc_discharge(self.wg['water_level'].value)
                            self.wg['water_level_out'].value = self.wg['water_level'].value
                        self.wg['wet_area'].value = cs.wet_area
                        self.wg['wet_perimeter'].value = cs.wet_perimeter
                        self.wg['hydraulic_radius'].value = cs.hydraulic_radius
                        self.wg['froude'].value = cs.froude_number
                        self.wg['freeboard_left'].value = cs.freeboard_left
                        self.wg['freeboard_right'].value = cs.freeboard_right
                        self._update_fig_geometry(cs)
                        if self.wg['water_level'].value <= min(yz_table.iloc[0, 1],yz_table.iloc[-1, 1]):
                            self.wg['warning'].value = ""
                        elif self.wg['water_level'].value <= max(yz_table.iloc[0, 1],yz_table.iloc[-1, 1]):
                            self.wg['warning'].value = "<p style='font-family:verdana; font-size:16px'>Profile is open on 1 side. The open side is taken as vertical (for area calculation). Perimeter and Hydraulic radius are based on wet profile (red line in graph). Adjust the profile table to correct</p>"
                        else:
                            self.wg['warning'].value = "<p style='font-family:verdana; font-size:16px'>Water level is higher than top of the profile. The sides above the profile are taken as vertical (for area calculation). Perimeter and Hydraulic radius are based on wet profile (red line in graph). Adjust the profile table to correct</p>"
                                 
                    else:
                        self.wg['warning'].value ="<p style='font-family:verdana; font-size:14px'>Profile provided should be non-decreasing. Please check file</p>"
                        self.wg['yz_message'].value ="<p style='font-family:verdana; font-size:18px'>Profile provided should be non-decreasing. Please check file</p>"
                        self.wg['discharge_out'].value = 0
                        self.wg['water_level_out'].value = 0
                        self.wg['wet_area'].value = 0
                        self.wg['wet_perimeter'].value = 0
                        self.wg['hydraulic_radius'].value = 0
                        self.wg['froude'].value = 0
                        self.wg['freeboard_left'].value = 0
                        self.wg['freeboard_right'].value = 0
        if self.wg['froude'].value == 0:
            self.wg['froude_message'].value = "<p style='font-family:verdana; font-size:14px'>No flow</p>"
        elif self.wg['froude'].value < 1:
            self.wg['froude_message'].value = "<p style='font-family:verdana; font-size:14px'>Flow is sub-critical</p>"
        elif self.wg['froude'].value > 1:
            self.wg['froude_message'].value = "<p style='font-family:verdana; font-size:14px'>Flow is super-critical</p>"
        else:
            self.wg['froude_message'].value = "<p style='font-family:verdana; font-size:14px'>Flow is critical</p>"        
    
    def _update_fig_geometry(self, cs):
        with self.figdwp.batch_update():
            if self.wg['cs_type'].value == 'Trapezoidal':
                self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['channel'])].x = cs._cs_x_coords()
                self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['channel'])].y = cs._cs_y_coords()
                
                self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['wet_profile'])].x = cs._cswet_x_coords()
                self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['wet_profile'])].y = cs._cswet_y_coords()
             
             # Water level
                self.figdwp.data[
                    get_pos_trace(self.figdwp, self.figure_geometry_legend['water level'])].x = [min(cs._cswet_x_coords()),
                                                                                             max(cs._cswet_x_coords())]
                self.figdwp.data[
                    get_pos_trace(self.figdwp, self.figure_geometry_legend['water level'])].y = [
                    self.wg['water_level_out'].value, self.wg['water_level_out'].value]
                self.figdwp.layout.xaxis.dtick = 10 ** np.ceil(np.log10(max(cs._cswet_x_coords())/10) )
           
            elif self.wg['cs_type'].value == 'YZ':
                self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['channel'])].x = cs.yz_values[:,0]
                self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['channel'])].y = cs.yz_values[:,1]
            
                self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['wet_profile'])].x = cs.wet_profile[:,0]
                self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['wet_profile'])].y = cs.wet_profile[:,1]

            # Water level
                self.figdwp.data[
                    get_pos_trace(self.figdwp, self.figure_geometry_legend['water level'])].x = [min(cs.wet_profile[:,0]),
                                                                                             max(cs.wet_profile[:,0])]
                self.figdwp.data[
                    get_pos_trace(self.figdwp, self.figure_geometry_legend['water level'])].y = [
                    self.wg['water_level_out'].value, self.wg['water_level_out'].value]
