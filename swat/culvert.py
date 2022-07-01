from math import acos, sqrt
import numpy as np
from math import cos, sin, pi
import warnings

class Gelok:
    """
    Basic formula of Gelok (1969) of hydraulic calculations of culverts, bridges and siphons
    Source: Cultuur technisch vademecum 1988 deel V par 1.1.6 (blz 803)
    """

    def __init__(self):
        pass

    def gelok_calculate_q(self, dz: float, wet_area: float, discharge_coeff: float, g: float = 9.81):
        """GELOK formula to calculate discharge
        Source with explanation: Cultuur technisch vademecum 1988 deel V par 1.1.6 (blz 803)

        :param dz: waterlevel difference over culvert [m]
        :param wet_area: wet area of the culver cross section [m2]
        :param discharge_coeff: discharge coefficient mu [-]
        :param g: gravitational acceleration [m/s^2]
        :return: discharge [m3/s]
        """
        return float(discharge_coeff * wet_area * sqrt(2 * g * dz))

    def gelok_calculate_dz(self, q: float, wet_area: float, discharge_coeff: float, g: float = 9.81):
        """GELOK formula to calculate water level difference
        Source with explanation: Cultuur technisch vademecum 1988 deel V par 1.1.6 (blz 803)

        :param q: discharge [m3/s]
        :param wet_area: wet area of the culver cross section [m2]
        :param discharge_coeff: discharge coefficient mu [-]
        :param g: gravitational acceleration [m/s^2]
        :return:  waterlevel difference over culvert [m]
        """
        return float(pow(q, 2) * pow(wet_area, -2)) / (pow(discharge_coeff, 2) * 2 * g)

    def _calculate_discharge_coeff(self, loss_coeff_in, loss_coeff_out, loss_coeff_friction, loss_coeff_additional):
        """
        Calculate friction coefficient from Loss coefficients
        :param loss_coeff_in: inlet loss coefficient, depends on inlet shape of the culvert [-]
        :param loss_coeff_out: outlet loss coefficient [-]
        :param loss_coeff_friction: friction loss coefficient
        :param loss_coeff_additional: additional loss coefficients
        :return:
        """
        return float(1 / sqrt(loss_coeff_in + loss_coeff_friction + loss_coeff_out + loss_coeff_additional))

    def _calculate_loss_coeff_friction(self, length:float, hydraulic_radius:float, k_manning:float, g: float=9.81):
        """
        Calculates loss coefficient due to friction
        """
        chezy = k_manning * pow(hydraulic_radius, 1 / 6)
        return float(2 * g * length / (pow(chezy, 2) * hydraulic_radius))

    def _calculate_loss_coeff_out(self, frac_wet_area: float, k_loss_outflow: float = 1):
        """Outlet loss coefficient (uittreeverlies)
        frac_wet_area: fraction wet area culvert of wet area channel
        k_loss_outflow: energy loss dependent from the outlet shape of the culvert (k=1: all energy is los, k=0: no energy loss)
        """
        return float(pow((1 - frac_wet_area), 2) * k_loss_outflow)

    def _hydraulic_radius(self, wet_perimeter: float, wet_area: float):
        """
        Calculates hydraulic radius
        """
        if wet_perimeter <= 0:
            warnings.warn(f"wet_perimeter is 0 or negative ({wet_perimeter}), returning None")
            return None
        if wet_area <= 0:
            warnings.warn(f"wet_area is 0 or negative ({wet_area}), returning None")
            return None
        return float(wet_area / wet_perimeter)

    def _velocity(self, q, wet_area):
        return q/wet_area

    def froude_number(self):
        raise NotImplementedError

    @staticmethod
    def xy_param_profile(param_profile, bob, height):
        # {'bottomwidth': 1, 'talud_l': 2, 'talud_r': 2}
        bottomwidth = param_profile['bottomwidth']
        talud_l_width = param_profile['talud_l'] * height
        talud_r_width  = param_profile['talud_r'] * height

        x = [-0.5 * bottomwidth - talud_l_width, -0.5 * bottomwidth,
             0.5 * bottomwidth, 0.5 * bottomwidth + talud_r_width]
        y = [height + bob, bob, bob, height + bob]
        return [x, y]



class CulvertRound(Gelok):
    def __init__(self, diameter: float, length: float, h_water: float, bob: float = 0, h_sediment: float = -999,
                 k_manning: float = 75, inlet_loss: float = 0.6, frac_wet_area_outflow: float = 0,
                 k_loss_outflow: float = 1, additional_loss: float = 0, param_profile: dict = None, g: float = 9.81):
        """
             Gelok formula for round culvert. Water level difference over culvert unknown.
             Source with explanation: Cultuur technisch vademecum 1988 deel V par 1.1.6 (blz 803)

             :param diameter: Culvert diameter [m]
             :param length: Culvert length [m]
             :param h_water: Absolute water level downstream of culvert [m+ref]
             :param bob: Invert level [m+ref]
             :param h_sediment: Absolute sediment level [m+ref]
             :param k_manning: K_manning (= 1/n) friction coefficient of culvert material [m^1/3 /s]
             :param inlet_loss: inlet loss coefficient, depends on inlet shape of the culvert [-]
             :param frac_wet_area_outflow: fraction of the wet area of the culvert compared to the wet area of the channel [-]
             :param k_loss_outflow: energy loss dependent from the outlet shape of the culvert (k=1: all energy is los, k=0: no energy loss)
             :param param_profile: parameterized downstream profile with the following keys: {'bottomwidth': 1, 'talud_l': 2, 'talud_r': 2}
             using this parameter will overwrite frac_wet_area_outflow
             :param additional_loss: additional energy losses due to e.g. corners [-]
             :param g: gravitational acceleration [m/s^2]

             """

        self.profile = {"diameter": diameter, "shape": "round"}
        # self.diameter = diameter
        self.length = length
        self.h_water_ds = h_water
        self.bob = bob
        self.h_sediment = h_sediment
        self.k_manning = k_manning
        self.cf_inlet_loss = inlet_loss
        self.frac_wet_area_outflow = frac_wet_area_outflow
        self.k_loss_outflow = k_loss_outflow
        self.cf_additional_loss = additional_loss
        self.param_profile = param_profile
        self.gravity = g


    def calc_dz(self, q: float):
        """
        Gelok formula for round culvert. Water level difference over culvert unknown.
        Source with explanation: Cultuur technisch vademecum 1988 deel V par 1.1.6 (blz 803)
        :param q: discharge [m3/s]
        :return: waterlevel difference over culvert [m]
        """
        if self.discharge_coefficient is not None:
            return self.gelok_calculate_dz(q, self.wet_area, self.discharge_coefficient, self.gravity)
        else:
            return None

    def calc_q(self, dz: float):
        """
        Gelok formula for round culvert. Discharge unknown.
        Source with explanation: Cultuur technisch vademecum 1988 deel V par 1.1.6 (blz 803)
        :param dz: waterlevel difference over culvert
        :return: discharge [m3/s]
        """
        if self.discharge_coefficient is not None:
            return self.gelok_calculate_q(dz, self.wet_area, self.discharge_coefficient, self.gravity)
        else:
            return None

    @property
    def wet_area(self):
        return self._wet_area_round(self.profile['diameter'], self.h_water_ds - self.bob, self.h_sediment - self.bob)

    @property
    def wet_perimeter(self):
        return self._wet_perimeter_round(self.profile['diameter'], self.h_water_ds - self.bob,
                                         self.h_sediment - self.bob)

    @property
    def hydraulic_radius(self):
        return self._hydraulic_radius(self.wet_perimeter, self.wet_area)

    @property
    def cf_friction_loss(self):
        if self.hydraulic_radius is not None:
            return self._calculate_loss_coeff_friction(self.length, self.hydraulic_radius, self.k_manning)
        else:
            return None

    @property
    def cf_outlet_loss(self):
        if self.param_profile is not None:
            water_depth = self.h_water_ds - self.bob
            width_upper = self.param_profile['bottomwidth'] + \
                          self.param_profile['talud_l'] * water_depth + \
                          self.param_profile['talud_r'] * water_depth
            profile_A = 0.5 * (self.param_profile['bottomwidth'] + width_upper) * water_depth
            # print('frac_wet_area_outflow has been re-determined using the parameterized profile input and the wet area of the culvert! '
            #       'If this is not intentional, please ensure to initiate the culvert with "param_profile = None"')
            self.frac_wet_area_outflow = self.wet_area / profile_A
        return self._calculate_loss_coeff_out(self.frac_wet_area_outflow, self.k_loss_outflow)

    @property
    def discharge_coefficient(self):
        if self.cf_friction_loss is not None:
            return self._calculate_discharge_coeff(self.cf_inlet_loss, self.cf_friction_loss, self.cf_outlet_loss,
                                                   self.cf_additional_loss)
        else:
            return None

    @property
    def xy_coords(self):
        return self._xy_coords_round(self.profile['diameter'], self.bob)


    def _wet_area_round(self, diameter: float, wl_above_bob: float, d_sediment: float=0):
        r'''Calculates the partial area of a circle
        .. math::
            \text{SA} = R^2\cos^{-1}\frac{(R - h)}{R} - (R - h)\sqrt{(2Rh - h^2)}

        To handle a layer of sediment as well as air, this function will calculate the area of sediment and the area of
        air, and will substract these areas from the full area of the circle. So:
        Wet area = circle area - sediment area - air area

        Parameters
        ----------
        D : float
            Diameter of the circle, [m]
        wl_above_bob : float
            Height measured from bottom of circle to liquid level, [m]
        d_sediment : float
            Height measured from bottom of circle to sediment level, [m]
        Returns
        -------
        SA_partial : float
            Partial (wetted) surface area, [m^2]
        Notes
        -----
        This method is undefined for :math:`h > D` and :math:`h < 0`, but those
        cases are handled by returning the full surface area and the zero
        respectively.

        References
        ----------
        .. [1] Weisstein, Eric W. "Circular Segment." Text. Wolfram Research, Inc.
           Accessed May 10, 2020. https://mathworld.wolfram.com/CircularSegment.html.
        '''
        if wl_above_bob > diameter:
            wl_above_bob = diameter  # Catch the case of a computed `h` being trivially larger than `D` due to floating point
        elif wl_above_bob < 0.0:
            return 0.0

        if d_sediment >= wl_above_bob:
            return 0
        if d_sediment < 0.0:
            d_sediment = 0

        R = 0.5 * diameter

        air_d = diameter - wl_above_bob

        area_air = R * R * acos((R - air_d) / R) - (R - air_d) * sqrt(2.0 * R * air_d - air_d * air_d)
        area_sediment = R * R * acos((R - d_sediment) / R) - (R - d_sediment) * \
                        sqrt(2.0 * R * d_sediment - d_sediment * d_sediment)

        area_whole_culvert = np.pi * R * R

        SA = area_whole_culvert - area_air - area_sediment
        if SA < 0.0:
            SA = 0.0  # Catch trig errors
        return SA

    def _wet_perimeter_round(self, diameter: float, wl_above_bob: float, d_sediment: float=0):
        r'''Calculates the wet perimeter of a circle.

        Parameters
        ----------
        D : float
            Diameter of the circle, [m]
        wl_above_bob : float
            Height measured from bottom of circle to liquid level, [m]
        d_sediment : float
            Height measured from bottom of circle to sediment level [m]
        Returns
        -------
        SA_partial : float
            Partial (wetted) surface area, [m^2]
        Notes
        -----
        This method is undefined for :math:`h > D` and :math:`h < 0`, but those
        cases are handled by returning the full surface area and the zero
        respectively.

        References
        ----------
        .. [1] Weisstein, Eric W. "Circular Segment." Text. Wolfram Research, Inc.
           Accessed May 10, 2020. https://mathworld.wolfram.com/CircularSegment.html.
        '''
        if wl_above_bob > diameter:
            wl_above_bob = diameter  # Catch the case of a computed `h` being trivially larger than `D` due to floating point
        if wl_above_bob < 0.0:
            return 0.0

        R = 0.5 * diameter
        P = 2 * R * acos((R - wl_above_bob) / R)

        if d_sediment > 0:
            # substract the arc length of the diameter submerged by sediment and add the width of the sediment to P
            p_sediment = 2 * R * acos((R - d_sediment) / R)
            width_sediment = sqrt(R ** 2 - (R - d_sediment) ** 2)
            P = P - p_sediment + width_sediment

        if P < 0.0:
            P = 0
        return P

    @classmethod
    def _xy_coords_round(cls, diameter, bob):
        steps = np.linspace(0, 2 * np.pi, 100)
        y = 0.5 * diameter * np.sin(steps) + bob + 0.5 * diameter
        x = 0.5 * diameter * np.cos(steps)
        return [x, y]


class CulvertRectangle(Gelok):
    def __init__(self, width: float, height: float, length: float, h_water: float, bob: float = 0,
                 h_sediment: float = -999, k_manning: float = 75, inlet_loss: float = 0.6,
                 frac_wet_area_outflow: float = 0, k_loss_outflow: float = 1, additional_loss: float = 0,
                 param_profile: dict = None, g: float = 9.81):
        """
        Gelok formula for round culvert. Water level difference over culvert unknown.
        Source with explanation: Cultuur technisch vademecum 1988 deel V par 1.1.6 (blz 803)
        :param diameter: Culvert diameter [m]
        :param length: Culvert length [m]
        :param h_water: Absolute water level downstream of culvert [m+ref]
        :param bob: Invert level [m+ref]
        :param h_sediment: Absolute water sediment level [m+ref]
        :param k_manning: K_manning (= 1/n) friction coefficient of culvert material [m^1/3 /s]
        :param inlet_loss: inlet loss coefficient, depends on inlet shape of the culvert [-]
        :param frac_wet_area_outflow: fraction of the wet area of the culvert compared to the wet area of the channel [-]
        :param k_loss_outflow: energy loss dependent from the outlet shape of the culvert (k=1: all energy is los, k=0: no energy loss)
        :param additional_loss: additional energy losses due to e.g. corners [-]
        :param param_profile: parameterized downstream profile with the following keys: {'bottomwidth': 1, 'talud_l': 2, 'talud_r': 2}
        using this parameter will overwrite frac_wet_area_outflow
        :param g: gravitational acceleration [m/s^2]

        """
        self.profile = {"width": width, "height": height, "shape": "rectangle"}
        # self.width = width
        # self.height = height
        self.length = length
        self.h_water_ds = h_water
        self.bob = bob
        self.h_sediment = h_sediment
        self.k_manning = k_manning
        self.cf_inlet_loss = inlet_loss
        self.frac_wet_area_outflow = frac_wet_area_outflow
        self.k_loss_outflow = k_loss_outflow
        self.cf_additional_loss = additional_loss
        self.param_profile = param_profile
        self.gravity = g

    def calc_dz(self, q: float):
        """
        Gelok formula for round culvert. Water level difference over culvert unknown.
        Source with explanation: Cultuur technisch vademecum 1988 deel V par 1.1.6 (blz 803)
        :param q: discharge [m3/s]
        :return: waterlevel difference over culvert [m]
        """
        if self.discharge_coefficient is not None:
            return self.gelok_calculate_dz(q, self.wet_area, self.discharge_coefficient, self.gravity)
        else:
            return None

    def calc_q(self, dz: float):
        """
        Gelok formula for round culvert. Discharge unknown.
        Source with explanation: Cultuur technisch vademecum 1988 deel V par 1.1.6 (blz 803)
        :param dz: waterlevel difference over culvert
        :return: discharge [m3/s]
        """
        if self.discharge_coefficient is not None:
            return self.gelok_calculate_q(dz, self.wet_area, self.discharge_coefficient, self.gravity)
        else:
            return None

    @property
    def wet_area(self):
        return self._wet_area_rectangle(self.profile['width'], self.profile['height'], self.h_water_ds - self.bob,
                                        self.h_sediment - self.bob)

    @property
    def wet_perimeter(self):
        return self._wet_perimeter_rectangle(self.profile['width'], self.profile['height'], self.h_water_ds - self.bob,
                                             self.h_sediment - self.bob)

    @property
    def hydraulic_radius(self):
        return self._hydraulic_radius(self.wet_perimeter, self.wet_area)

    @property
    def cf_friction_loss(self):
        if self.hydraulic_radius is not None:
            return self._calculate_loss_coeff_friction(self.length, self.hydraulic_radius, self.k_manning)
        else:
            return None

    @property
    def cf_outlet_loss(self):
        if self.param_profile is not None:
            water_depth = self.h_water_ds - self.bob
            width_upper = self.param_profile['bottomwidth'] + \
                          self.param_profile['talud_l'] * water_depth + \
                          self.param_profile['talud_r'] * water_depth
            profile_A = 0.5 * (self.param_profile['bottomwidth'] + width_upper) * water_depth
            # print('frac_wet_area_outflow has been re-determined using the parameterized profile input and the wet area of the culvert! '
            #       'If this is not intentional, please ensure to initiate the culvert with "param_profile = None"')
            self.frac_wet_area_outflow = self.wet_area / profile_A
        return self._calculate_loss_coeff_out(self.frac_wet_area_outflow, self.k_loss_outflow)

    @property
    def discharge_coefficient(self):
        if self.cf_friction_loss is not None:
            return self._calculate_discharge_coeff(self.cf_inlet_loss, self.cf_friction_loss, self.cf_outlet_loss,
                                                   self.cf_additional_loss)
        else:
            return None

    @property
    def xy_coords(self):
        return self._xy_coords_rectangle(self.profile['width'], self.profile['height'], self.bob)


    def _wet_area_rectangle(self, width: float, height: float, wl_above_bob: float, d_sediment: float = 0):
        r'''Calculates the wet area of a rectangle.
        .. math::
            \text{SA} = R^2\cos^{-1}\frac{(R - h)}{R} - (R - h)\sqrt{(2Rh - h^2)}
        Parameters
        ----------
        D : float
            Diameter of the circle, [m]
        wl_above_bob : float
            Height measured from bottom of rectangle to liquid level, [m]
        d_sediment : float
            Height measured from bottom of rectangle to liquid level, [m]
        Returns
        -------
        SA_partial : float
            Partial (wetted) surface area, [m^2]
        Notes
        -----
        This method is undefined for :math:`h > D` and :math:`h < 0`, but those
        cases are handled by returning the full surface area and the zero
        respectively.

        References
        ----------
        .. [1] Weisstein, Eric W. "Circular Segment." Text. Wolfram Research, Inc.
           Accessed May 10, 2020. https://mathworld.wolfram.com/CircularSegment.html.
        '''
        if wl_above_bob > height:
            wl_above_bob = height  # Catch the case of a computed `h` being trivially larger than `D` due to floating point
        if wl_above_bob < 0.0:
            return 0.

        if d_sediment >= wl_above_bob:
            return 0
        if d_sediment < 0.0:
            d_sediment = 0

        SA = (wl_above_bob - d_sediment) * width
        if SA < 0.0:
            SA = 0.0  # Catch trig errors
        return SA


    def _wet_perimeter_rectangle(self, width:float, height: float, wl_above_bob: float, d_sediment: float=0):
        r'''Calculates the wet perimeter of a rectangle.

        Parameters
        ----------
        D : float
            Diameter of the circle, [m]
        wl_above_bob : float
            Height measured from bottom of circle to liquid level, [m]
        d_sediment : float
            Height measured from bottom of rectangle to liquid level, [m]
        Returns
        -------
        P : float
            Wetted perimeter, [m]
        '''
        if wl_above_bob < 0.0:
            return 0.0

        if d_sediment >= wl_above_bob:
            return 0
        if d_sediment < 0.0:
            d_sediment = 0

        if wl_above_bob > height:
            P = (height - d_sediment) * 2 + width * 2
        else:
            P = (wl_above_bob - d_sediment) * 2 + width

        if P < 0.0:
            P = 0
        return P


    @classmethod
    def _xy_coords_rectangle(cls, width:float, heigth:float, bob):
        x = [-0.5*width, 0.5*width, 0.5*width, -0.5*width, -0.5*width]
        y = [bob, bob, bob+heigth, bob+heigth, bob]

        return [x, y]


if __name__ == "__main__":
    width = 0.50
    cv = CulvertRectangle(width=width, height=1.2, length=10, h_water=-0.3, bob=-0.2, h_sediment=-999,
                                         k_manning=75, inlet_loss=0.6, frac_wet_area_outflow=0,
                                         k_loss_outflow=1, additional_loss=0)
    print(cv.calc_dz(0.2))


