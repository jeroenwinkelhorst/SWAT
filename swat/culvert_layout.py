from swat.culvert import CulvertRound, CulvertRectangle, Gelok
from swat.channel import Trapezoidal
import ipywidgets as widgets
from ipywidgets.widgets import HBox, VBox
import plotly.graph_objs as go
import numpy as np
import copy


def get_pos_trace(fig, legend_name):
    """Return position of trace in fig.data based on legend name"""
    for i in range(len(fig.data)):
        if fig.data[i].name == legend_name:
            return i


class CulvertNotebookLayout:
    def __init__(self):
        self.output_widget = widgets.Output()
        self.figure_widget = widgets.Output()
        self.input_widgets = {}
        self.output_widgets = {}
        self.figure_widgets = {}
        self.figure_geometry_legend = {}
        self.widget_style = {}

        self._set_general_vars()

    def __repr__(self):
        return "Culvert notebook layout"

    def _set_general_vars(self):
        self.widget_style['fixed_width'] = {'description_width': '200px'}

    def run(self):
        self._initial_field_definitions()
        return VBox([
            widgets.HTML(value="<p style='font-family:verdana; font-size:24px; padding:20px'>Initial parameters</p>"),
            HBox([self.input_widgets['input_unknown'], widgets.HTML(
                value="<p style='font-family:verdana; font-size:12px;'>(dz = headloss; q = discharge)</p>")]),
            self.input_widgets['input_shape'],
            self.input_widgets['input_button1'],
            HBox([self.output_widget,self.figure_widget])
        ])

    def _warning_incorrect_result(self, warn: bool = True):
        if warn == True:
            self.output_widgets[
                'output_warning'].value = "<p style='font-family:verdana; font-size:24px; padding:20px; color:red'>Error in results - check your input!</p>"
        else:
            self.output_widgets['output_warning'].value = ""

    def _main_layout(self, b):
        self._main_field_definitions()
        self._create_fig_geometry()

        self.input_widgets['input_unknown'].disabled = True
        self.input_widgets['input_shape'].disabled = True

        # Hydraulic input
        with self.output_widget:
            display(
                widgets.HTML(
                    value="<p style='font-family:verdana; font-size:24px;padding:20px'>Hydraulic parameters</p>"))
            display(self.input_widgets['input_dz']) if self.input_widgets['input_unknown'].value == 'q' else None
            display(self.input_widgets['input_q']) if self.input_widgets['input_unknown'].value == 'dz' else None
            display(self.input_widgets['input_h'])

        # Selection Shape specific geometry
        # Input parameters based on Culvert shape
        if self.input_widgets['input_shape'].value == 'round':
            geometry_input = self._input_fields_culvert_round()
        elif self.input_widgets['input_shape'].value == 'rectangle':
            geometry_input = self._input_fields_culvert_rectangle()
        else:
            with output:
                display(widgets.Label("No culvert shape selected. Please run the cell again."))

        #  Geometry in/output
        with self.output_widget:
            display(
                VBox([
                    widgets.HTML(value="<p style='font-family:verdana; font-size:24px;padding:20px'>Geometry</p>"),
                    HBox([
                        VBox(geometry_input + [
                            self.input_widgets['input_length'],
                            self.input_widgets['input_bob'],
                        ]),
                        VBox([self.output_widgets['output_wet_area'],
                              self.output_widgets['output_wet_area_channel']
                              ]),
                    ]),
                    HBox([self.input_widgets['input_k_manning'], widgets.HTML(
                        value="<p style='font-family:verdana; font-size:12px;'> Roughness coefficient applies to overall profile.</p>")]),
                    HBox([self.input_widgets['input_h_sediment'], widgets.HTML(
                        value="<p style='font-family:verdana; font-size:12px;'> Consider adjusting roughness if preforming calculation with much sediment.</p>")]),

                    HBox([self.input_widgets['input_checkbox_profile']]),
                    HBox([self.input_widgets['input_profile_talud_l']]),
                    HBox([self.input_widgets['input_profile_talud_r']]),
                    HBox([self.input_widgets['input_profile_bottomwidth']]),
                    self.output_widgets['output_warning_profile'],
                    
                    widgets.HTML(value="<p style='font-family:verdana; font-size:24px;padding:20px'>Energy losses</p>"),
                    widgets.HTML(
                        value="<p style='font-family:verdana; font-size:12px;'>The parameter description can be found in de appendix</p>"),
                    VBox([
                        widgets.HTML(value="<p style='font-family:verdana; font-style: italic;'>Inlet loss</p>"),
                        self.input_widgets['input_inlet_loss'],
                        widgets.HTML(value="<p style='font-family:verdana; font-style: italic;'>Friction loss</p>"),
                        HBox([self.input_widgets['input_k_manning'], self.output_widgets['output_friction_loss']]),
                        widgets.HTML(value="<p style='font-family:verdana; font-style: italic;'>Outlet loss</p>"),
                        HBox([VBox([self.input_widgets['input_frac_wet_area_outflow'],
                                    self.input_widgets['input_k_loss_outflow'], ]),
                              self.output_widgets['output_outlet_loss']]),
                        widgets.HTML(value="<p style='font-family:verdana; font-style: italic;'>Additional loss</p>"),
                        self.input_widgets['input_additional_loss'],
                        widgets.HTML(
                            value="<p style='font-family:verdana; font-style: italic;'>Discharge coeffient</p>"),
                        self.output_widgets['output_discharge_coefficient']
                    ]),
                ]),
            )
        with self.figure_widget:
            display(self.figure_widgets['geometry'])
        # Results
        output_par = None
        output_par = self.output_widgets['output_q'] if self.input_widgets['input_unknown'].value == 'q' else output_par
        output_par = self.output_widgets['output_dz'] if self.input_widgets[
                                                             'input_unknown'].value == 'dz' else output_par

        with self.output_widget:
            display(widgets.HTML(value="<p style='font-family:verdana; font-size:24px;padding:20px'>Results</p>"))
            if self.input_widgets['input_unknown'].value == 'q':
                display(self.output_widgets['output_unknown_q'])
            elif self.input_widgets['input_unknown'].value == 'dz':
                display(self.output_widgets['output_unknown_dz'])
            display(self.output_widgets['output_warning'])
            display(self.input_widgets['input_sensitivity_add'])

    def _input_value_change(self, b):
        self._warning_incorrect_result(False)
        # Input parameters based on Unknown
        if self.input_widgets['input_shape'].value == 'round':
            if self.input_widgets['input_checkbox_profile'].value == False:
                dk = CulvertRound(
                    diameter=self.input_widgets['input_diameter'].value,
                    length=self.input_widgets['input_length'].value,
                    h_water=self.input_widgets['input_h'].value,
                    bob=self.input_widgets['input_bob'].value,
                    h_sediment=self.input_widgets['input_h_sediment'].value,
                    k_manning=self.input_widgets['input_k_manning'].value,
                    inlet_loss=self.input_widgets['input_inlet_loss'].value,
                    frac_wet_area_outflow=self.input_widgets['input_frac_wet_area_outflow'].value,
                    k_loss_outflow=self.input_widgets['input_k_loss_outflow'].value,
                    additional_loss=self.input_widgets['input_additional_loss'].value,
                    param_profile=None)
            else:
                dk = CulvertRound(
                    diameter=self.input_widgets['input_diameter'].value,
                    length=self.input_widgets['input_length'].value,
                    h_water=self.input_widgets['input_h'].value,
                    bob=self.input_widgets['input_bob'].value,
                    h_sediment=self.input_widgets['input_h_sediment'].value,
                    k_manning=self.input_widgets['input_k_manning'].value,
                    inlet_loss=self.input_widgets['input_inlet_loss'].value,
                    frac_wet_area_outflow=self.input_widgets['input_frac_wet_area_outflow'].value,
                    k_loss_outflow=self.input_widgets['input_k_loss_outflow'].value,
                    additional_loss=self.input_widgets['input_additional_loss'].value,
                    param_profile={'bottomwidth': self.input_widgets['input_profile_bottomwidth'].value,
                                   'talud_l': self.input_widgets['input_profile_talud_l'].value,
                                   'talud_r': self.input_widgets['input_profile_talud_r'].value})
            if self.input_widgets['input_unknown'].value == 'q':
                try:
                    self.output_widgets['output_q'].value = dk.calc_q(self.input_widgets['input_dz'].value)
                except:
                    self.output_widgets['output_q'].value = -999
                    self._warning_incorrect_result(True)

            elif self.input_widgets['input_unknown'].value == 'dz':
                try:
                    self.output_widgets['output_dz'].value = dk.calc_dz(self.input_widgets['input_q'].value)
                except:
                    self.output_widgets['output_dz'].value = -999
                    self._warning_incorrect_result(True)

        elif self.input_widgets['input_shape'].value == 'rectangle':
            if self.input_widgets['input_checkbox_profile'].value == False:
                dk = CulvertRectangle(
                    width=self.input_widgets['input_width'].value,
                    height=self.input_widgets['input_height'].value,
                    length=self.input_widgets['input_length'].value,
                    h_water=self.input_widgets['input_h'].value,
                    bob=self.input_widgets['input_bob'].value,
                    h_sediment=self.input_widgets['input_h_sediment'].value,
                    k_manning=self.input_widgets['input_k_manning'].value,
                    inlet_loss=self.input_widgets['input_inlet_loss'].value,
                    frac_wet_area_outflow=self.input_widgets['input_frac_wet_area_outflow'].value,
                    k_loss_outflow=self.input_widgets['input_k_loss_outflow'].value,
                    additional_loss=self.input_widgets['input_additional_loss'].value,
                    param_profile=None)
            else:
                dk = CulvertRectangle(
                    width=self.input_widgets['input_width'].value,
                    height=self.input_widgets['input_height'].value,
                    length=self.input_widgets['input_length'].value,
                    h_water=self.input_widgets['input_h'].value,
                    bob=self.input_widgets['input_bob'].value,
                    h_sediment=self.input_widgets['input_h_sediment'].value,
                    k_manning=self.input_widgets['input_k_manning'].value,
                    inlet_loss=self.input_widgets['input_inlet_loss'].value,
                    frac_wet_area_outflow=self.input_widgets['input_frac_wet_area_outflow'].value,
                    k_loss_outflow=self.input_widgets['input_k_loss_outflow'].value,
                    additional_loss=self.input_widgets['input_additional_loss'].value,
                    param_profile={'bottomwidth': self.input_widgets['input_profile_bottomwidth'].value,
                                   'talud_l': self.input_widgets['input_profile_talud_l'].value,
                                   'talud_r': self.input_widgets['input_profile_talud_r'].value})
            if self.input_widgets['input_unknown'].value == 'q':
                try:
                    self.output_widgets['output_q'].value = dk.calc_q(self.input_widgets['input_dz'].value)
                except:
                    self.output_widgets['output_q'].value = -999
                    self._warning_incorrect_result(True)

            elif self.input_widgets['input_unknown'].value == 'dz':
                try:
                    self.output_widgets['output_dz'].value = dk.calc_dz(self.input_widgets['input_q'].value)
                except:
                    self.output_widgets['output_dz'].value = -999
                    self._warning_incorrect_result(True)

        try:
            self.output_widgets['output_wet_area'].value = dk.wet_area
            # General output values
            self.output_widgets['output_discharge_coefficient'].value = dk.discharge_coefficient
            self.output_widgets['output_friction_loss'].value = dk.cf_friction_loss
            self.output_widgets['output_outlet_loss'].value = dk.cf_outlet_loss
        except:
            self.output_widgets['output_wet_area'].value = -999
            self.output_widgets['output_discharge_coefficient'].value = -999
            self.output_widgets['output_friction_loss'].value = -999
            self.output_widgets['output_outlet_loss'].value = -999
            self._warning_incorrect_result(True)

        self.output_widgets[
            'output_unknown_q'].value = "<p style='font-family:verdana; font-size:14px;'>Given the head loss of " + str(
            self.input_widgets['input_dz'].value) + "m, the resulting discharge is " + str(
            round(self.output_widgets['output_q'].value, 3)) + " m3/s</p>"
        self.output_widgets[
            'output_unknown_dz'].value = "<p style='font-family:verdana; font-size:14px;'>Given the discharge of " + str(
            self.input_widgets['input_q'].value) + "m3/s, the resulting headloss is " + str(
            round(self.output_widgets['output_dz'].value, 3)) + " m</p>"

        self._update_fig_geometry(dk)

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

        # Culvert
        legend_culvert = "culvert"
        fig.add_scatter(x=CulvertRound._xy_coords_round(1, 0)[0],
                        y=CulvertRound._xy_coords_round(1, 0)[1], line=dict(color="black"),
                        name=legend_culvert)
        # Water level
        legend_water_level = "water level downstream"
        x_min = CulvertRound._xy_coords_round(1, 0)[0].min()
        x_max = CulvertRound._xy_coords_round(1, 0)[0].max()
        fig.add_scatter(x=[x_min, x_max], y=[self.input_widgets['input_h'].value, self.input_widgets['input_h'].value],
                        line=dict(color="blue"),
                        name=legend_water_level)

        # Sediment level
        legend_sediment_level = "Sediment level"
        x_min = CulvertRound._xy_coords_round(1, 0)[0].min()
        x_max = CulvertRound._xy_coords_round(1, 0)[0].max()
        fig.add_scatter(x=[x_min, x_max], y=[self.input_widgets['input_h_sediment'].value,
                                             self.input_widgets['input_h_sediment'].value],
                        line=dict(color="brown"),
                        name=legend_sediment_level)

        # channel
        legend_channel = "channel"

        self.figure_geometry_legend['culvert'] = legend_culvert
        self.figure_geometry_legend['water level'] = legend_water_level
        self.figure_geometry_legend['sediment level'] = legend_sediment_level
        self.figure_geometry_legend['channel'] = legend_channel
        self.figure_widgets['geometry'] = fig

    def _update_fig_geometry(self, dk):
        with self.figure_widgets['geometry'].batch_update():
            x_min = min(dk.xy_coords[0])
            x_max = max(dk.xy_coords[0])
            # Profile
            if self.input_widgets['input_checkbox_profile'].value == True:
                profile_x_vals, profile_y_vals = Gelok.xy_param_profile(
                    param_profile={'bottomwidth': self.input_widgets['input_profile_bottomwidth'].value,
                                   'talud_l': self.input_widgets['input_profile_talud_l'].value,
                                   'talud_r': self.input_widgets['input_profile_talud_r'].value},
                    bob=self.input_widgets['input_bob'].value, height=self.input_widgets['input_h'].value)
                # profile_x_vals = [float(i) for i in self.input_widgets['input_profile_x'].value.strip().split()]
                # profile_y_vals = [float(i) for i in self.input_widgets['input_profile_y'].value.strip().split()]
                if (len(profile_x_vals) == len(profile_y_vals)) and len(profile_x_vals) >= 2:
                    self.output_widgets['output_warning_profile'].value = ""
                    self.figure_widgets['geometry'].data[get_pos_trace(self.figure_widgets['geometry'],
                                                                       self.figure_geometry_legend[
                                                                           'channel'])].x = profile_x_vals
                    self.figure_widgets['geometry'].data[get_pos_trace(self.figure_widgets['geometry'],
                                                                       self.figure_geometry_legend[
                                                                           'channel'])].y = profile_y_vals
                    # wet_x_coords = list(zip(
                    #     *calculate_wet_profile(profile_x_vals, profile_y_vals, self.input_widgets['input_h'].value)))[0]
                    x_min = profile_x_vals[0]
                    x_max = profile_x_vals[-1]
                    self.output_widgets['output_wet_area_channel'].value = Trapezoidal(
                        k_manning=20, bottom_slope=0.001,
                        slope_left = self.input_widgets['input_profile_talud_l'].value,
                        slope_right=self.input_widgets['input_profile_talud_r'].value,
                        bed_width=self.input_widgets['input_profile_bottomwidth'].value,
                        bed_level=self.input_widgets['input_bob'].value,
                        surface_level=self.input_widgets['input_h'].value)._wet_area(self.input_widgets['input_h'].value)
                    # self.output_widgets['output_wet_area_channel'].value = calculate_wet_area(profile_x_vals,
                    #                                                                           profile_y_vals,
                    #                                                                           self.input_widgets[
                    #                                                                               'input_h'].value)
                    self.input_widgets['input_frac_wet_area_outflow'].value = self.output_widgets[
                                                                                  'output_wet_area'].value / \
                                                                              self.output_widgets[
                                                                                  'output_wet_area_channel'].value
                else:
                    self.output_widgets[
                        'output_warning_profile'].value = "ERROR: X and Y coordinates are invalid - please correct your input"

                    # Water level
            self.figure_widgets['geometry'].data[
                get_pos_trace(self.figure_widgets['geometry'], self.figure_geometry_legend['water level'])].x = [
                x_min, x_max]
            self.figure_widgets['geometry'].data[
                get_pos_trace(self.figure_widgets['geometry'], self.figure_geometry_legend['water level'])].y = [
                self.input_widgets['input_h'].value, self.input_widgets['input_h'].value]

                    # Sediment level
            self.figure_widgets['geometry'].data[
                get_pos_trace(self.figure_widgets['geometry'], self.figure_geometry_legend['sediment level'])].x = [
                x_min, x_max]
            self.figure_widgets['geometry'].data[
                get_pos_trace(self.figure_widgets['geometry'], self.figure_geometry_legend['sediment level'])].y = [
                self.input_widgets['input_h_sediment'].value, self.input_widgets['input_h_sediment'].value]

            # Culvert
            self.figure_widgets['geometry'].data[
                get_pos_trace(self.figure_widgets['geometry'], self.figure_geometry_legend['culvert'])].x = \
                dk.xy_coords[0]
            self.figure_widgets['geometry'].data[
                get_pos_trace(self.figure_widgets['geometry'], self.figure_geometry_legend['culvert'])].y = \
                dk.xy_coords[1]


    def _input_fields_culvert_round(self):
        style = self.widget_style['fixed_width']
        self.input_widgets['input_diameter'] = widgets.BoundedFloatText(
            value=1,
            min=0.01,
            max=10.0,
            step=0.01,
            description='Culvert diameter [m]:',
            disabled=False,
            style=style
        )
        self.input_widgets['input_diameter'].observe(self._input_value_change, names='value')
        return [self.input_widgets['input_diameter']]

    def _input_fields_culvert_rectangle(self):
        style = self.widget_style['fixed_width']
        self.input_widgets['input_width'] = widgets.BoundedFloatText(
            value=1,
            min=0.01,
            max=10.0,
            step=0.01,
            description='Culvert width [m]:',
            disabled=False,
            style=style
        )

        self.input_widgets['input_height'] = widgets.BoundedFloatText(
            value=1,
            min=0.01,
            max=10.0,
            step=0.01,
            description='Culvert height [m]:',
            disabled=False,
            style=style
        )
        self.input_widgets['input_width'].observe(self._input_value_change, names='value')
        self.input_widgets['input_height'].observe(self._input_value_change, names='value')
        return [self.input_widgets['input_width'], self.input_widgets['input_height']]

    def _input_checkbox_profile_change(self, b):
        self.input_widgets['input_profile_talud_l'].disabled = False if self.input_widgets[
                                                                      'input_checkbox_profile'].value == True else True
        self.input_widgets['input_profile_talud_r'].disabled = False if self.input_widgets[
                                                                      'input_checkbox_profile'].value == True else True
        self.input_widgets['input_profile_bottomwidth'].disabled = False if self.input_widgets[
                                                                      'input_checkbox_profile'].value == True else True
        # self.input_widgets['input_profile_x'].disabled = False if self.input_widgets[
        #                                                               'input_checkbox_profile'].value == True else True

        if self.input_widgets['input_checkbox_profile'].value == True:
            self.input_widgets['input_frac_wet_area_outflow'].disabled = True
            self.input_widgets['input_checkbox_profile'].disabled = True
            self.figure_widgets['geometry'].add_scatter(x=[-1, 0, 1], y=[1, 0, 1], line=dict(color="green"),
                                                        name=self.figure_geometry_legend['channel'])
            self._input_value_change(1)

    def _add_sensitivity_plot(self, click):
        self._sensitivity_field_definitions()
        self.figure_widgets['sensitivity'] = go.FigureWidget()
        self.input_widgets['input_sensitivity_add'].disabled = True

        self.figure_widgets['sensitivity'].layout.template = "plotly_white"
        self.figure_widgets['sensitivity'].add_scatter(line=dict(color="black"))
        self._update_sensitivity_plot(1)

        with self.output_widget:
            display(VBox([
                self.input_widgets['input_sensitivity'],
                self.figure_widgets['sensitivity'],
                self.input_widgets['input_sensitivity_axis'],
            ])
            )

    def _update_sensitivity_plot_slider(self, change):
        self._update_sensitivity_plot(1, slider_input=True)

    def _update_sensitivity_plot(self, change, slider_input=False):

        # Define x-axis parameter
        if self.input_widgets['input_sensitivity'].value == 'dz':
            self.figure_widgets['sensitivity'].update_xaxes(title_text='Head Loss [m]')
            sens_x = self.input_widgets['input_dz']
        elif self.input_widgets['input_sensitivity'].value == 'q':
            self.figure_widgets['sensitivity'].update_xaxes(title_text='Discharge [m3/s]')
            sens_x = self.input_widgets['input_q']
        elif self.input_widgets['input_sensitivity'].value == 'k_manning':
            self.figure_widgets['sensitivity'].update_xaxes(title_text='K Manning [m^1/3 /s]')
            sens_x = self.input_widgets['input_k_manning']
        elif self.input_widgets['input_sensitivity'].value == 'diameter':
            self.figure_widgets['sensitivity'].update_xaxes(title_text='Diameter [m]')
            sens_x = self.input_widgets['input_diameter']
        elif self.input_widgets['input_sensitivity'].value == 'width':
            self.figure_widgets['sensitivity'].update_xaxes(title_text='Width [m]')
            sens_x = self.input_widgets['input_width']
        elif self.input_widgets['input_sensitivity'].value == 'height':
            self.figure_widgets['sensitivity'].update_xaxes(title_text='Height [m]')
            sens_x = self.input_widgets['input_height']
        elif self.input_widgets['input_sensitivity'].value == 'length':
            self.figure_widgets['sensitivity'].update_xaxes(title_text='Length [m]')
            sens_x = self.input_widgets['input_length']
        elif self.input_widgets['input_sensitivity'].value == 'water level':
            self.figure_widgets['sensitivity'].update_xaxes(title_text='Water level [m+ref]')
            sens_x = self.input_widgets['input_h']
        else:
            raise ValueError('Input parameter unknown not defined')

        # Define range x-axis
        if slider_input:
            min_sens_x = self.input_widgets['input_sensitivity_axis'].value[0]
            max_sens_x = self.input_widgets['input_sensitivity_axis'].value[1]
        else:
            n = 50
            min_sens_x = sens_x.value - n * sens_x.step if sens_x.value - n * sens_x.step > sens_x.min else sens_x.min
            max_sens_x = sens_x.value + n * sens_x.step if sens_x.value + n * sens_x.step < sens_x.max else sens_x.max
            self.input_widgets['input_sensitivity_axis'].value = [min_sens_x, max_sens_x]
            self.input_widgets['input_sensitivity_axis'].min = sens_x.min
            self.input_widgets['input_sensitivity_axis'].max = sens_x.max
            self.input_widgets['input_sensitivity_axis'].step = sens_x.step

        def get_y_for_xrange_unknown_q(obj, property_to_vary, x_list):
            """
            obj: Culvert object
            property_to_vary: parameter to vary
            x_list: coordinates of the x-axis
            """
            tmp_obj = copy.deepcopy(obj)
            result = []
            for i in x_list:
                if property_to_vary in tmp_obj.profile.keys():
                    tmp_obj.profile[property_to_vary] = i
                else:
                    setattr(tmp_obj, property_to_vary, i)
                result.append(tmp_obj.calc_q(self.input_widgets['input_dz'].value))
            return result

        def get_y_for_xrange_unknown_dz(obj, property_to_vary, x_list):
            """
            obj: Culvert object
            property_to_vary: parameter to vary
            x_list: coordinates of the x-axis
            """
            tmp_obj = copy.deepcopy(obj)
            result = []
            for i in x_list:
                if property_to_vary in tmp_obj.profile.keys():
                    tmp_obj.profile[property_to_vary] = i
                else:
                    setattr(tmp_obj, property_to_vary, i)
                result.append(tmp_obj.calc_dz(self.input_widgets['input_q'].value))
            return result

        # Culvert definition
        if self.input_widgets['input_shape'].value == 'round':
            if self.input_widgets['input_checkbox_profile'].value == False:
                dk_sens = CulvertRound(
                    diameter=self.input_widgets['input_diameter'].value,
                    length=self.input_widgets['input_length'].value,
                    h_water=self.input_widgets['input_h'].value,
                    bob=self.input_widgets['input_bob'].value,
                    h_sediment=self.input_widgets['input_h_sediment'].value,
                    k_manning=self.input_widgets['input_k_manning'].value,
                    inlet_loss=self.input_widgets['input_inlet_loss'].value,
                    frac_wet_area_outflow=self.input_widgets[
                        'input_frac_wet_area_outflow'].value,
                    k_loss_outflow=self.input_widgets[
                        'input_k_loss_outflow'].value,
                    additional_loss=self.input_widgets[
                        'input_additional_loss'].value)
            else:
                dk_sens = CulvertRound(
                    diameter=self.input_widgets['input_diameter'].value,
                    length=self.input_widgets['input_length'].value,
                    h_water=self.input_widgets['input_h'].value,
                    bob=self.input_widgets['input_bob'].value,
                    h_sediment=self.input_widgets['input_h_sediment'].value,
                    k_manning=self.input_widgets['input_k_manning'].value,
                    inlet_loss=self.input_widgets['input_inlet_loss'].value,
                    frac_wet_area_outflow=self.input_widgets['input_frac_wet_area_outflow'].value,
                    k_loss_outflow=self.input_widgets['input_k_loss_outflow'].value,
                    additional_loss=self.input_widgets['input_additional_loss'].value,
                    param_profile={'bottomwidth': self.input_widgets['input_profile_bottomwidth'].value,
                                   'talud_l': self.input_widgets['input_profile_talud_l'].value,
                                   'talud_r': self.input_widgets['input_profile_talud_r'].value})
        elif self.input_widgets['input_shape'].value == 'rectangle':
            if self.input_widgets['input_checkbox_profile'].value == False:
                dk_sens = CulvertRectangle(
                    width=self.input_widgets['input_width'].value,
                    height=self.input_widgets['input_height'].value,
                    length=self.input_widgets['input_length'].value,
                    h_water=self.input_widgets['input_h'].value,
                    bob=self.input_widgets['input_bob'].value,
                    h_sediment=self.input_widgets['input_h_sediment'].value,
                    k_manning=self.input_widgets['input_k_manning'].value,
                    inlet_loss=self.input_widgets['input_inlet_loss'].value,
                    frac_wet_area_outflow=self.input_widgets[
                        'input_frac_wet_area_outflow'].value,
                    k_loss_outflow=self.input_widgets[
                        'input_k_loss_outflow'].value,
                    additional_loss=self.input_widgets[
                        'input_additional_loss'].value)
            else:
                dk_sens = CulvertRectangle(
                    width=self.input_widgets['input_width'].value,
                    height=self.input_widgets['input_height'].value,
                    length=self.input_widgets['input_length'].value,
                    h_water=self.input_widgets['input_h'].value,
                    bob=self.input_widgets['input_bob'].value,
                    h_sediment=self.input_widgets['input_h_sediment'].value,
                    k_manning=self.input_widgets['input_k_manning'].value,
                    inlet_loss=self.input_widgets['input_inlet_loss'].value,
                    frac_wet_area_outflow=self.input_widgets['input_frac_wet_area_outflow'].value,
                    k_loss_outflow=self.input_widgets['input_k_loss_outflow'].value,
                    additional_loss=self.input_widgets['input_additional_loss'].value,
                    param_profile={'bottomwidth': self.input_widgets['input_profile_bottomwidth'].value,
                                   'talud_l': self.input_widgets['input_profile_talud_l'].value,
                                   'talud_r': self.input_widgets['input_profile_talud_r'].value})
        else:
            raise ValueError('Invalid shape definition')

        # Define y-values
        sens_x_axis = np.linspace(min_sens_x, max_sens_x, 100)
        sens_y_axis = []

        if self.input_widgets['input_unknown'].value == 'dz':
            self.figure_widgets['sensitivity'].update_yaxes(title_text='Head Loss [m]')
            if self.input_widgets['input_sensitivity'].value == 'dz':
                sens_y_axis = sens_x_axis
            elif self.input_widgets['input_sensitivity'].value == 'q':
                sens_y_axis = [dk_sens.calc_dz(i) for i in sens_x_axis]
            else:
                sens_y_axis = get_y_for_xrange_unknown_dz(dk_sens, self.input_widgets['input_sensitivity'].value,
                                                          sens_x_axis)
        elif self.input_widgets['input_unknown'].value == 'q':
            self.figure_widgets['sensitivity'].update_yaxes(title_text='Discharge [m3/s]')
            if self.input_widgets['input_sensitivity'].value == 'q':
                sens_y_axis = sens_x_axis
            elif self.input_widgets['input_sensitivity'].value == 'dz':
                sens_y_axis = [dk_sens.calc_q(i) for i in sens_x_axis]
            else:
                sens_y_axis = get_y_for_xrange_unknown_q(dk_sens, self.input_widgets['input_sensitivity'].value,
                                                         sens_x_axis)
        else:
            raise ValueError("input_unknown incorrect", input_unknown.value)

        self.figure_widgets['sensitivity'].data[0].x = sens_x_axis
        self.figure_widgets['sensitivity'].data[0].y = sens_y_axis

    def _initial_field_definitions(self):
        self.input_widgets['input_unknown'] = widgets.ToggleButtons(
            options=['dz', 'q'],
            description='Parameter to calculate:',
            disabled=False,
            button_style='',  # 'success', 'info', 'warning', 'danger' or ''
            tooltips=['Head loss', 'Discharge'],
            style=self.widget_style['fixed_width']
        )

        self.input_widgets['input_shape'] = widgets.ToggleButtons(
            options=['round', 'rectangle'],
            description='Choose culvert shape:',
            disabled=False,
            button_style='',  # 'success', 'info', 'warning', 'danger' or ''
            tooltips=[''],
            style=self.widget_style['fixed_width']
        )

        self.input_widgets['input_button1'] = widgets.Button(
            description='Freeze initial parameters and start calculation',
            disabled=False,
            button_style='',  # 'success', 'info', 'warning', 'danger' or ''
            tooltip='',
            icon='',
            layout=widgets.Layout(width='initial', margin='20px')
        )
        self.input_widgets['input_button1'].on_click(self._main_layout)

    def _main_field_definitions(self):
        style = self.widget_style['fixed_width']
        self.input_widgets['input_q'] = widgets.BoundedFloatText(
            value=1,
            min=0,
            max=10.0,
            step=0.01,
            description='Discharge [m3/s]:',
            disabled=False,
            style=style
        )

        self.input_widgets['input_dz'] = widgets.BoundedFloatText(
            value=0.02,
            min=0,
            max=1,
            step=0.001,
            description='Head loss [m]:',
            disabled=False,
            style=style
        )

        self.input_widgets['input_length'] = widgets.BoundedFloatText(
            value=10,
            min=0,
            max=1000.0,
            step=0.1,
            description='Culvert length [m]:',
            disabled=False,
            style=style
        )

        self.input_widgets['input_h'] = widgets.BoundedFloatText(
            value=0.75,
            min=-1000,
            max=10000.0,
            step=0.01,
            description='Water level [m+ref]:',
            disabled=False,
            style=style
        )

        self.input_widgets['input_bob'] = widgets.BoundedFloatText(
            value=0,
            min=-1000,
            max=10000.0,
            step=0.01,
            description='Invert level [m+ref]:',
            disabled=False,
            style=style
        )

        self.input_widgets['input_h_sediment'] = widgets.BoundedFloatText(
            value=0,
            min=-1000,
            max=10000.0,
            step=0.01,
            description='Sediment level [m+ref]:',
            disabled=False,
            style=style
        )

        self.input_widgets['input_k_manning'] = widgets.BoundedIntText(
            value=75,
            min=1,
            max=200,
            step=1,
            description=r'Roughness \(k_M [ m^{1/3} / s ]\):',
            disabled=False,
            style=style
        )

        self.input_widgets['input_inlet_loss'] = widgets.BoundedFloatText(
            value=0.6,
            min=0,
            max=1,
            step=0.01,
            description=r'Inlet loss  \(\varepsilon_{in} \) [-]:',
            disabled=False,
            style=style
        )

        self.input_widgets['input_frac_wet_area_outflow'] = widgets.BoundedFloatText(
            value=0.5,
            min=0,
            max=1,
            step=0.01,
            description=r'Fraction wet area \(\alpha \ [-] \):',
            disabled=False,
            style=style
        )

        self.input_widgets['input_k_loss_outflow'] = widgets.BoundedFloatText(
            value=1,
            min=0,
            max=1,
            step=0.01,
            description=r'Shape coefficient \( k [-] \):',
            disabled=False,
            style=style
        )

        self.input_widgets['input_additional_loss'] = widgets.BoundedFloatText(
            value=0,
            min=0,
            max=10.0,
            step=0.01,
            description=r'\(\varepsilon_{other} \) [-]:',
            disabled=False,
            style=style
        )

        self.input_widgets['input_checkbox_profile'] = widgets.Checkbox(
            value=False,
            description='Add channel profile (downstream)',
            disabled=False,
            indent=False,
            style=style
        )
        self.input_widgets['input_checkbox_profile'].observe(self._input_checkbox_profile_change, names='value')

        self.input_widgets['input_profile_talud_l'] = widgets.BoundedFloatText(
            value=1,
            min=0.01,
            max=100,
            step=0.01,
            description='Channel - Left slope [m/m]',
            disabled=True,
            style=style
        )

        self.input_widgets['input_profile_talud_r'] = widgets.BoundedFloatText(
            value=1,
            min=0.01,
            max=100,
            step=0.01,
            description='Channel - Right slope [m/m]',
            disabled=True,
            style=style
        )

        self.input_widgets['input_profile_bottomwidth'] = widgets.BoundedFloatText(
            value=0,
            min=0,
            max=100,
            step=0.01,
            description='Channel - Bottom width [m]',
            disabled=True,
            style=style
        )

        # self.input_widgets['input_profile_y'] = widgets.Text(
        #     value='1 0 1',
        #     placeholder='space separated',
        #     description='Channel profile y-coords:',
        #     disabled=True,
        #     style=style
        # )

        self.input_widgets['input_sensitivity_add'] = widgets.Button(
            description='Add sensitivity plot',
            disabled=False,
            button_style='',  # 'success', 'info', 'warning', 'danger' or ''
            tooltip='',
            icon='',
            style=style,
            layout=widgets.Layout(margin='20px'))

        self.output_widgets['output_q'] = widgets.FloatText(description="Discharge [m3/s]", style=style, disabled=True)
        self.output_widgets['output_dz'] = widgets.FloatText(description="Water level difference [m]", style=style,
                                                             disabled=True)
        self.output_widgets['output_wet_area'] = widgets.FloatText(description="Wet Area culvert [m2]", style=style,
                                                                   disabled=True)
        self.output_widgets['output_wet_area_channel'] = widgets.FloatText(description="Wet Area channel [m2]",
                                                                           style=style, disabled=True)
        self.output_widgets['output_discharge_coefficient'] = widgets.FloatText(
            description=r"Discharge coefficient \( \mu \) [-]",
            style=style,
            disabled=True)
        self.output_widgets['output_friction_loss'] = widgets.FloatText(
            description=r"Friction loss \( \varepsilon_{friction} \) [-]",
            style=style, disabled=True)
        self.output_widgets['output_outlet_loss'] = widgets.FloatText(
            description=r"Outlet loss \( \varepsilon_{out} \) [-]",
            style=style, disabled=True)

        self.output_widgets['output_froude'] = widgets.FloatText(description="Froude number []", style=style,
                                                                 disabled=True)
        self.output_widgets['output_warning'] = widgets.HTML(value=f"", style=style)
        self.output_widgets['output_warning_profile'] = widgets.Label(value="", style=style)

        self.output_widgets['output_unknown_q'] = widgets.HTML(
            value="<p style='font-family:verdana; font-size:14px;'>Given the head loss of " + str(
                self.input_widgets['input_dz'].value) + "m, the resulting discharge is " + str(
                round(self.output_widgets['output_q'].value, 3)) + " m3/s</p>")
        self.output_widgets['output_unknown_dz'] = widgets.HTML(
            value="<p style='font-family:verdana; font-size:14px;'>Given the discharge of " + str(
                self.input_widgets['input_q'].value) + "m3/s, the resulting headloss is " + str(
                round(self.output_widgets['output_dz'].value, 3)) + " m</p>")

        self.input_widgets['input_q'].observe(self._input_value_change, names='value')
        self.input_widgets['input_dz'].observe(self._input_value_change, names='value')
        self.input_widgets['input_length'].observe(self._input_value_change, names='value')
        self.input_widgets['input_h'].observe(self._input_value_change, names='value')
        self.input_widgets['input_bob'].observe(self._input_value_change, names='value')
        self.input_widgets['input_h_sediment'].observe(self._input_value_change, names='value')
        self.input_widgets['input_k_manning'].observe(self._input_value_change, names='value')
        self.input_widgets['input_inlet_loss'].observe(self._input_value_change, names='value')
        self.input_widgets['input_frac_wet_area_outflow'].observe(self._input_value_change, names='value')
        self.input_widgets['input_k_loss_outflow'].observe(self._input_value_change, names='value')
        self.input_widgets['input_additional_loss'].observe(self._input_value_change, names='value')
        self.input_widgets['input_profile_talud_l'].observe(self._input_value_change, names='value')
        self.input_widgets['input_profile_talud_r'].observe(self._input_value_change, names='value')
        self.input_widgets['input_profile_bottomwidth'].observe(self._input_value_change, names='value')
        # self.input_widgets['input_profile_x'].observe(self._input_value_change, names='value')
        # self.input_widgets['input_profile_y'].observe(self._input_value_change, names='value')
        self.input_widgets['input_sensitivity_add'].on_click(self._add_sensitivity_plot)

    def _sensitivity_field_definitions(self):
        input_par_list = ['q', 'dz', 'k_manning']
        if self.input_widgets['input_shape'].value == 'round':
            input_par_list.extend(['diameter'])
        if self.input_widgets['input_shape'].value == 'rectangle':
            input_par_list.extend(['width', 'height'])

        self.input_widgets['input_sensitivity'] = widgets.ToggleButtons(
            options=input_par_list,
            description='X-axis value:',
            disabled=False,
            button_style='',  # 'success', 'info', 'warning', 'danger' or ''
            tooltips=['Discharge', 'Water level difference over culvert'],
            #             style=style
        )

        self.input_widgets['input_sensitivity_axis'] = widgets.FloatRangeSlider(
            value=[0.5, 1],
            min=0.1,
            max=10.0,
            step=0.1,
            description='Range x-axis:',
            disabled=False,
            continuous_update=False,
            orientation='horizontal',
            readout=True,
            readout_format='.3f',
            style={'description_width': '80px'}
        )
        self.input_widgets['input_sensitivity_axis'].layout.width = '800px'

        self.input_widgets['input_sensitivity'].observe(self._update_sensitivity_plot)
        self.input_widgets['input_sensitivity_axis'].observe(self._update_sensitivity_plot_slider)


test = CulvertNotebookLayout()
test.run()

if __name__ == "__main__":
    test = CulvertNotebookLayout()
    test.run()
