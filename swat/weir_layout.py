import ipywidgets as widgets
from ipywidgets.widgets import HBox, VBox
import plotly.graph_objs as go
from swat.weir import Weir
import numpy as np

def get_pos_trace(fig, legend_name):
    """Return position of trace in fig.data based on legend name"""
    for i in range(len(fig.data)):
        if fig.data[i].name == legend_name:
            return i

class WeirNotebookLayout:
    def __init__(self):
        self.output = widgets.Output()
        self.figure = widgets.Output()
        # widgets dict
        self.wg = {}
        self.wg_style = {}

        self.figure_geometry_legend = {}

        self._set_general_vars()
        self._widget_defs()
        self._create_fig_geometry()

    def __repr__(self):
        return "Weir notebook layout"

    def _set_general_vars(self):
        self.wg_style['fixed_width'] = {'description_width': '400px'}
        self.wg_style['layout'] = {'width': '550px'}


    def run(self):
        return VBox([
            widgets.HTML(value="<p style='font-family:verdana; font-size:24px; padding:20px'>Initial parameters</p>"),
            self.wg['weir_type'],
            HBox([self.output,self.figure])])

    def _add_general_fields(self):
        with self.output:
            display(
                VBox([
                    self.wg['water_level_up'],
                    self.wg['water_level_down']
                ])
            )
        with self.figure:
            display(self.figdwp)

    def _add_weir_fields(self):
        with self.output:
            display(
                VBox([
                    widgets.HTML(
                        value="<p style='font-family:verdana; font-size:24px; padding:20px'>Weir geometry</p>"),
                    self.wg['crest_level'],
                    self.wg['width'],
                    self.wg['lateral_coef'],
                    self.wg['weir_coef'],
                    self.wg['loss_coef'],
                    
                ])
            )

    def _add_orifice_fields(self):
        with self.output:
            display(
                VBox([
                    widgets.HTML(
                        value="<p style='font-family:verdana; font-size:24px; padding:20px'>Orifice geometry</p>"),
                    self.wg['crest_level'],
                    self.wg['gate_level'],
                    self.wg['width'],
                    self.wg['lateral_coef'],
                    self.wg['weir_coef'],
                    self.wg['loss_coef'],
                    self.wg['contract_coef']
                ])
            )
                    

    def _add_result_fields(self):
        with self.output:
            display(
                VBox([
                    widgets.HTML(value="<p style='font-family:verdana; font-size:24px; padding:20px'>Output</p>"),
                    self.wg['flow_type'],
                    self.wg['discharge_out'],
                    self.wg['wet_area'],
                    self.wg['velocity'],
                    self.wg['hac'],
                ])
            )

    def _observe_weir_type(self, b):
        self.wg['weir_type'].disabled = True
        self._add_general_fields()
        if self.wg['weir_type'].value == "weir":
            self.wg['gate_level'].value = np.inf
            self._add_weir_fields()
        if self.wg['weir_type'].value == "orifice":
            self._add_orifice_fields()    
        self._add_result_fields()
        self._update_fig_geometry()

    def _observe_weir_input(self, b):
        self._update()

    def _widget_defs(self):
        style = self.wg_style['fixed_width']
        layout = self.wg_style['layout']

        self.wg['weir_type'] = widgets.Dropdown(
            options=['', 'weir','orifice'],
            value='',
            description='Select weir type:',
            disabled=False,
            style=style,
            layout = layout
        )

        self.wg['weir_type'].observe(self._observe_weir_type, names='value')

        # General input
        self.wg['water_level_up'] = widgets.BoundedFloatText(
            value=2,
            min=-10,
            max=200,
            step=0.01,
            description='upstream water level [m+datum]',
            disabled=False,
            style=style,
            layout = layout
        )
        self.wg['water_level_up'].observe(self._observe_weir_input, names='value')

        self.wg['water_level_down'] = widgets.BoundedFloatText(
            value=2,
            min=-10,
            max=200,
            step=0.01,
            description='downstream water level [m+datum]',
            disabled=False,
            style=style,
            layout = layout
        )
        self.wg['water_level_down'].observe(self._observe_weir_input, names='value')

        # ouput
        self.wg['flow_type'] = widgets.Dropdown(
            options=['', 'free','submerged'],
            value='',
            description='flow type:',
            disabled=True,
            style=style,
            layout = layout
        )
        self.wg['discharge_out'] = widgets.BoundedFloatText(
            value=0,
            min=-2000000,
            max=2000000,
            description='Discharge [m3/s]:',
            disabled=True,
            style=style,
            layout = layout
        )

        self.wg['wet_area'] = widgets.BoundedFloatText(
            value=0,
            max=2000000,            
            description='Wet area over crest [m2]:',
            disabled=True,
            style=style,
            layout = layout
        )
        self.wg['velocity'] = widgets.BoundedFloatText(
            value=0,
            min=-2000000,
            max=2000000,            
            description='Velocity over crest [m/s]:',
            disabled=True,
            style=style,
            layout = layout
        )       
        self.wg['hac'] = widgets.BoundedFloatText(
            value=1,
            min=0,
            max=2000,            
            description='Level above crest (overstortende straal) [m]:',
            disabled=True,
            style=style,
            layout = layout
        )
        # geometry simple
        self.wg['crest_level'] = widgets.BoundedFloatText(
            value=1,
            min=-10,
            max=2000,
            step=0.01,
            description='crest level [m+datum]:',
            disabled=False,
            style=style,
            layout = layout
        )
        self.wg['crest_level'].observe(self._observe_weir_input, names='value')

        self.wg['gate_level'] = widgets.BoundedFloatText(
            value=1,
            min=0,
            max=2000,
            step=0.01,
            description='gate level above crest [m]:',
            disabled=False,
            style=style,
            layout = layout
        )
        
        self.wg['gate_level'].observe(self._observe_weir_input, names='value')
      
        self.wg['width'] = widgets.BoundedFloatText(
            value=1,
            min=0,
            max=500,
            step=0.01,
            description='Width [m]:',
            disabled=False,
            style=style,
            layout = layout
        )
        
        self.wg['width'].observe(self._observe_weir_input, names='value')
        
        self.wg['lateral_coef'] = widgets.BoundedFloatText(
            value=1,
            min=0,
            max=1,
            step=0.01,
            description='lateral contraction coefficient [-]:',
            disabled=False,
            style=style,
            layout = layout
        )
        
        self.wg['lateral_coef'].observe(self._observe_weir_input, names='value')
        
        self.wg['contract_coef'] = widgets.BoundedFloatText(
            value=0.63,
            min=0,
            max=1,
            step=0.01,
            description='gate contraction coefficient (standard 0.63) [-]:',
            disabled=False,
            style=style,
            layout = layout
        )
        
        self.wg['contract_coef'].observe(self._observe_weir_input, names='value')        
        
        self.wg['weir_coef'] = widgets.BoundedFloatText(
            value=1,
            min=0,
            max=1.5,
            step=0.01,
            description='weir coeffcient (depends on crest shape, only used for free flow)',
            disabled=False,
            style=style,
            layout = layout
        )
        
        self.wg['weir_coef'].observe(self._observe_weir_input, names='value')         

        self.wg['loss_coef'] = widgets.BoundedFloatText(
            value=0.9,
            min=0,
            max=1,
            step=0.01,
            description='loss_coeffcient (used with submerged flow)',
            disabled=False,
            style=style,
            layout = layout
        )
        
        self.wg['loss_coef'].observe(self._observe_weir_input, names='value')
    
    def _create_fig_geometry(self):
        # Figure widget
        fig = go.FigureWidget()
        fig.layout.template = "plotly_white"
        #fig.layout.xaxis.dtick = 1
        fig.update_xaxes()
        fig.update_xaxes(title_text='[m]')
        fig.update_yaxes(title_text='[m]')
        fig.update_yaxes(range=[0.1, 2.5])  
        fig.update_xaxes(range=[0.5, 1.5])
        
        
        legend_wl_critical = "critical water level"
        fig.add_scatter(x=[1, 2], y=[2, 2],
                        line=dict(color="#000000", width=1 ),
                        name=legend_wl_critical, mode='lines')
        legend_wl_up = "water level upstream"
        fig.add_scatter(x=[0, 1], y=[2, 2],
                        line=dict(color="#000090", width=2 ),
                        name=legend_wl_up, mode='lines') 
        legend_wl_down = "water level downstream"
        fig.add_scatter(x=[1, 2], y=[2, 2],
                        line=dict(color="#6060FF", width=2),
                        name=legend_wl_down, mode='lines')
        # channel
        legend_weir = "Weir"
        fig.add_scatter(x=[1, 1], y=[0, 1],
                        line=dict(color="green", width=10),
                        name=legend_weir, mode='lines')

        legend_gate = "Gate"
        fig.add_scatter(x=[1, 1], y=[3, 2],
                        line=dict(color="red", width=10),
                        name=legend_gate, mode='lines')
                # Water level

        self.figure_geometry_legend['weir'] = legend_weir
        self.figure_geometry_legend['gate'] = legend_gate
        self.figure_geometry_legend['wl_up'] = legend_wl_up
        self.figure_geometry_legend['wl_down'] = legend_wl_down
        self.figure_geometry_legend['wl_critical'] = legend_wl_critical
        self.figdwp = fig

    def _update(self):
        self._update_fig_geometry()
        if self.wg['weir_type'].value == 'weir':
            
            struc = Weir(
                crest = self.wg['crest_level'].value,
                width = self.wg['width'].value,
                Cw = self.wg['weir_coef'].value,
                Cs = self.wg['loss_coef'].value,
                Cl = self.wg['lateral_coef'].value,
                water_level_up = self.wg['water_level_up'].value,
                water_level_down = self.wg['water_level_down'].value)
               
            if struc.flow_check:
                self.wg['flow_type'].value = 'free'
            else:
                self.wg['flow_type'].value = 'submerged'
            self.wg['discharge_out'].value = struc.discharge
            self.wg['wet_area'].value = struc.wet_area
            self.wg['velocity'].value = struc.velocity
            self.wg['hac'].value = struc.height_above_crest 

           # self._update_fig_geometry(cs)
        elif self.wg['weir_type'].value == 'orifice':
            
            struc = Weir(
                crest = self.wg['crest_level'].value,
                width = self.wg['width'].value,
                gate = self.wg['gate_level'].value,
                Cc = self.wg['contract_coef'].value,
                Cw = self.wg['weir_coef'].value,
                Cs = self.wg['loss_coef'].value,
                Cl = self.wg['lateral_coef'].value,
                water_level_up = self.wg['water_level_up'].value,
                water_level_down = self.wg['water_level_down'].value)      
            if struc.flow_check:
                self.wg['flow_type'].value = 'free'
            else:
                self.wg['flow_type'].value = 'submerged'
            self.wg['discharge_out'].value = struc.discharge
            self.wg['wet_area'].value = struc.wet_area
            self.wg['velocity'].value = struc.velocity
            self.wg['hac'].value = struc.height_above_crest

    def _update_fig_geometry(self):
        with self.figdwp.batch_update():
            if self.wg['weir_type'].value == 'orifice':
                minval = min( self.wg['crest_level'].value, self.wg['water_level_up'].value, self.wg['water_level_down'].value)
                maxval = max( self.wg['crest_level'].value + self.wg['gate_level'].value, self.wg['water_level_up'].value, self.wg['water_level_down'].value)
                deltarange = max(0.25 * (maxval - minval), 0.5)
                if self.wg['water_level_up'].value >= self.wg['water_level_down'].value:
                    critical = min(max(self.wg['water_level_up'].value - self.wg['crest_level'].value, 0) / 1.5, self.wg['gate_level'].value) + self.wg['crest_level'].value
                    self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['wl_critical'])].x = [1, 2]
                else:
                    critical = min(max(self.wg['water_level_down'].value-self.wg['crest_level'].value,0) / 1.5,self.wg['gate_level'].value) + self.wg['crest_level'].value
                    self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['wl_critical'])].x = [0, 1]    

                self.figdwp.update_yaxes(range=[minval-deltarange,maxval+deltarange])

                self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['weir'])].y = [minval-deltarange, self.wg['crest_level'].value]
                self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['gate'])].y = [self.wg['crest_level'].value + self.wg['gate_level'].value, maxval+deltarange]
                self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['wl_up'])].y = [self.wg['water_level_up'].value, self.wg['water_level_up'].value]
                self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['wl_down'])].y = [self.wg['water_level_down'].value, self.wg['water_level_down'].value]
                self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['wl_critical'])].y = [critical, critical]
              
           
            elif self.wg['weir_type'].value == 'weir':
                minval = min( self.wg['crest_level'].value, self.wg['water_level_up'].value, self.wg['water_level_down'].value)
                maxval = max( self.wg['crest_level'].value, self.wg['water_level_up'].value, self.wg['water_level_down'].value)
                deltarange = max(0.25 * (maxval - minval), 0.5)
                if self.wg['water_level_up'].value >= self.wg['water_level_down'].value:
                    critical = max(self.wg['water_level_up'].value - self.wg['crest_level'].value,0) / 1.5 + self.wg['crest_level'].value
                    self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['wl_critical'])].x = [1, 2]
                else:
                    critical = max(self.wg['water_level_down'].value - self.wg['crest_level'].value,0) / 1.5 + self.wg['crest_level'].value
                    self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['wl_critical'])].x = [0, 1] 

                self.figdwp.update_yaxes(range=[minval - deltarange,maxval + deltarange])

                self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['weir'])].y = [minval-deltarange, self.wg['crest_level'].value]
                self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['gate'])].y = [maxval+deltarange + 2, maxval + deltarange + 2]
                self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['wl_up'])].y = [self.wg['water_level_up'].value, self.wg['water_level_up'].value]
                self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['wl_down'])].y = [self.wg['water_level_down'].value, self.wg['water_level_down'].value]
                self.figdwp.data[get_pos_trace(self.figdwp, self.figure_geometry_legend['wl_critical'])].y = [critical, critical]

