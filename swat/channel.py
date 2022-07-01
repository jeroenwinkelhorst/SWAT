import numpy as np
import math

class Manning:

    def _calc_discharge(self, k_manning: float, slope: float, wet_area: float,
                        hydraulic_radius: float):
        return wet_area * hydraulic_radius ** (2 / 3) * k_manning * slope ** (1 / 2)

    def _velocity(self, discharge, wet_area):
        if wet_area>0:
            return discharge / wet_area
        else:
            return 0    

    def _chezy_to_k_manning(self, k_manning):
        pass

    def bos_en_bijkerk_to_k_manning(self, k_manning,depth):
        """
        Calculate friction value Km based on Bos en Bijkerk
        """
        return k_manning*depth**(1/3)
        

   

class Trapezoidal(Manning):
    def __init__(self, k_manning: float, bottom_slope: float, slope_left: float, slope_right: float, bed_width: float,
                 bed_level: float = 0, surface_level: float=1, b_b: bool=False):
        """
        Calculate discharge through a Trapizoidal channel.
        k_manning = friction value (ks strickler)
        bottom_slope = longitudinal slope of the channel
        slope_left = left side slope angle
        slope_right =  right side slope angle
        bed_width =  bed width
        bed_level =  bed level in m above datum
        surface level = surface level that determines top of the profile side slopes
        b_b = boolean that determines if Bos en Bijkerk should be used
        """
        self.slope_left = slope_left
        self.slope_right = slope_right
        self.bed_width = bed_width
        self.bed_level = bed_level
        self.k_manning = k_manning
        self.b_b= b_b
        self.bottom_slope = bottom_slope
        self.surface_level = surface_level

        # parameters dependent on input
        self.water_level = None
        self.discharge = None

    @property
    def wet_area(self):
        if self.water_level is not None:
            return self._wet_area(self.water_level)
        else:
            raise ValueError("Water level unknown, perform calculation")



    def _wet_area(self, water_level):
        depth = max(min(water_level, self.surface_level) - self.bed_level, 0)
        flood = max(water_level - max(self.surface_level, self.bed_level), 0)
        return self.bed_width * depth + 0.5 * depth * (depth * self.slope_left) + 0.5 * depth * (
                depth * self.slope_right) + (self.bed_width + depth * self.slope_left + depth * self.slope_right) * flood

    @property
    def depth(self):
        if self.water_level is not None:
            return self.water_level - self.bed_level
        else:
            raise ValueError("Water level unknown, perform calculation")

    @property
    def wet_perimeter(self):
        if self.water_level is not None:
            return self._calculate_wet_perimeter(self.water_level)
        else:
            raise ValueError("Water level unknown, perform calculation")

    @property
    def hydraulic_radius(self):
        if self.water_level is not None:
            return self._hydraulic_radius()
        else:
            raise ValueError("Water level unknown, perform calculation")

    @property
    def froude_number(self):
        if self.water_level is not None:
            return self._calc_froude()
        else:
            raise ValueError("Water level unknown, perform calculation")

    @property
    def freeboard(self):
        if self.water_level is not None:
            return self._calc_freeboard()
        else:
            raise ValueError("Water level unknown, perform calculation")

    @property
    def velocity(self):
        if self.water_level is not None:
            return self._velocity(self.discharge, self.wet_area)
        else:
            raise ValueError("Water level unknown, perform calculation")

    def _hydraulic_radius(self):
        return self.wet_area / self.wet_perimeter

    def _calculate_wet_perimeter(self, water_level: float):
        depth = max(min(water_level, self.surface_level) - self.bed_level,0)
        perimeter = (self.bed_width + math.sqrt(depth**2 + (depth * self.slope_left)**2) 
                     + math.sqrt( depth ** 2 + (depth * self.slope_right) ** 2))
        return perimeter        
    
    def _calc_froude(self):
        if self.wet_area>0:
            return self.discharge/(self.wet_area*math.sqrt(9.81*max(self.water_level - self.bed_level,0)))
        else:
            return 0

    def _calc_freeboard(self):
        return self.surface_level-self.water_level    

    def calc_discharge(self, water_level):
        self._reset_calculation_output()
        self.water_level = water_level
        if self.b_b:
            self.k_bb=self.bos_en_bijkerk_to_k_manning(self.k_manning,max(self.depth,0))
            self.discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
        else:
            self.discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius)
        return self.discharge

    def calc_water_level(self, discharge):
        """
        Iteratively determines which water level occurs at a given discharge. First, the discharge at a water depth of
        1 meter is determined. The water level is then either increased or decreased until the right discharge is found.
        """
        self._reset_calculation_output()
        self.discharge = discharge
        self.water_level = self.bed_level+1
        delta=0.1
        if self.b_b:
            self.k_bb=self.bos_en_bijkerk_to_k_manning(self.k_manning,max(self.depth,0))
            check_discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
        else:
            check_discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius)
        if check_discharge > discharge:
            delta = -delta
            while check_discharge > discharge:
                self.water_level += delta
                if self.b_b:
                    self.k_bb = self.bos_en_bijkerk_to_k_manning(self.k_manning,max(self.depth,0))
                    check_discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
                else:
                    check_discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            delta = -delta*0.1
            while check_discharge < discharge:
                self.water_level += delta
                if self.b_b:
                    self.k_bb = self.bos_en_bijkerk_to_k_manning(self.k_manning,max(self.depth,0))
                    check_discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
                else:
                    check_discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            delta = -delta * 0.1
            while check_discharge > discharge:
                self.water_level += delta
                if self.b_b:
                    self.k_bb = self.bos_en_bijkerk_to_k_manning(self.k_manning,max(self.depth,0))
                    check_discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
                else:
                    check_discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            self.water_level -= delta*0.5
            if self.b_b:
                self.k_bb=self.bos_en_bijkerk_to_k_manning(self.k_manning,max(self.depth,0))
                check_discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            else:
                check_discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius) 
            if check_discharge>discharge:
                self.water_level+=delta*0.5
            else:
                self.water_level-=delta*0.5
        elif check_discharge < discharge:
            while check_discharge < discharge:
                self.water_level += delta
                if self.b_b:
                    self.k_bb = self.bos_en_bijkerk_to_k_manning(self.k_manning,max(self.depth,0))
                    check_discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
                else:
                    check_discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            delta = -delta*0.1
            while check_discharge > discharge:
                self.water_level += delta
                if self.b_b:
                    self.k_bb = self.bos_en_bijkerk_to_k_manning(self.k_manning,max(self.depth,0))
                    check_discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
                else:
                    check_discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            delta = -delta*0.1
            while check_discharge < discharge:
                self.water_level += delta
                if self.b_b:
                    self.k_bb = self.bos_en_bijkerk_to_k_manning(self.k_manning,max(self.depth,0))
                    check_discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
                else:
                    check_discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            self.water_level -= delta*0.5
            if self.b_b:
                self.k_bb = self.bos_en_bijkerk_to_k_manning(self.k_manning,max(self.depth,0))
                check_discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            else:
                check_discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            if check_discharge < discharge:
                self.water_level += delta*0.5
            else:
                self.water_level -= delta*0.5                    
        return self.water_level 

    def _reset_calculation_output(self):
        self.water_level = None
        self.discharge = None

    def plot(self):
        pass

    def _cs_x_coords(self):
        if self.surface_level > self.bed_level:
            dx = [0, (self.surface_level - self.bed_level) * self.slope_left, self.bed_width,
              (self.surface_level - self.bed_level) * self.slope_right
              ]
        else:
            dx = [0, self.bed_width]  
        return list(np.cumsum(dx))

    def _cs_y_coords(self):
        if self.surface_level > self.bed_level:
            return [self.surface_level, self.bed_level, self.bed_level, self.surface_level]
        else:
             return [self.bed_level, self.bed_level]    

    def _cswet_x_coords(self):
        if self.surface_level > self.bed_level:
            dx = [min(max(0,(self.surface_level - self.water_level) * self.slope_left ), (self.surface_level - self.bed_level) * self.slope_left),
              (self.surface_level - self.bed_level) * self.slope_left,
              self.bed_width + (self.surface_level - self.bed_level) * self.slope_left,
              max((min(self.surface_level, self.water_level) - self.bed_level) * self.slope_right, 0) + self.bed_width + (self.surface_level - self.bed_level) * self.slope_left
              ]
        else:
            dx = [0,self.bed_width]
        return dx

    def _cswet_y_coords(self):
        if self.surface_level > self.bed_level:
            return [min(max(self.water_level, self.bed_level), self.surface_level), self.bed_level, self.bed_level, max(min(self.water_level, self.surface_level), self.bed_level)]
        else:
             return [ self.bed_level, self.bed_level]

class YZ(Manning):
    def __init__(self, k_manning: float, bottom_slope: float, yz_values: np.array, b_b: bool=False):
        
        self.k_manning = k_manning
        self.b_b= b_b
        self.bottom_slope = bottom_slope
        self.yz_values = yz_values

        # parameters saved by functions
        self.water_level = None
        self.discharge = None

    def yz_ok(self):
        ok=True
        for i in range(len(self.yz_values)-1):
            if self.yz_values[i, 0] > self.yz_values[i + 1, 0]:
                ok=False
                continue
        return ok

    @property
    def wet_area(self):
        if self.water_level is not None:
            if self.water_level <= min(self.yz_values[:,1]):
                return 0
            else:    
                return self._calculate_wet_area(self.water_level)
        else:
            raise ValueError("Water level unknown, perform calculation")

    @property
    def hydraulic_radius(self):
        if self.water_level is not None:
            if self.water_level <= min(self.yz_values[:,1]):
                return 0
            else:
                return self.wet_area/self.wet_perimeter
        else:
            raise ValueError("Water level unknown, perform calculation")

    @property
    def wet_profile(self):
        if self.water_level is not None:
            if self.water_level <= min(self.yz_values[:,1]):
                return np.array([[0,self.water_level]])
            else:
                return self._calculate_wet_profile(self.water_level)
        else:
            raise ValueError("Water level unknown, perform calculation")

    @property
    def wet_perimeter(self):
        if self.water_level is not None:
            return self._calculate_wet_perimeter(self.water_level)
        else:
            raise ValueError("Water level unknown, perform calculation")

    @property
    def froude_number(self):
        if self.water_level is not None:
            return self._calc_froude()
        else:
            raise ValueError("Water level unknown, perform calculation")

    @property
    def freeboard_left(self):
        if self.water_level is not None:
            return self._calc_freeboard_left()
        else:
            raise ValueError("Water level unknown, perform calculation")

    @property
    def freeboard_right(self):
        if self.water_level is not None:
            return self._calc_freeboard_right()
        else:
            raise ValueError("Water level unknown, perform calculation")

    def _calculate_wet_area(self, h: float):
        return abs(np.trapz(x=self.wet_profile[:, 0], y=self.wet_profile[:, 1] - h))

    def _calculate_wet_profile(self, h: float):
        """Calculate profile (coordinates of a channel (Trapezoidal method)"""
    
        yz_coords_new = []
        for i in range(len(self.yz_values)):
            if self.yz_values[i,1] > h and self.yz_values[i-1,1] <= h and i > 0:
                yz_coords_new.append(
                    (
                        (h - self.yz_values[i - 1,1])
                        / (self.yz_values[i, 1] - self.yz_values[i - 1 ,1])
                        * (self.yz_values[i, 0] - self.yz_values[i - 1 ,0])
                        + self.yz_values[i - 1,0],
                        h,
                    )
                )
            elif self.yz_values[i,1] <= h:
                yz_coords_new.append(self.yz_values[i,])    
            if i != len(self.yz_values) - 1:
                if self.yz_values[i,1] > h and self.yz_values[i + 1, 1] <= h:
                    yz_coords_new.append(
                        (
                            (h - self.yz_values[i, 1])
                            / (self.yz_values[i + 1, 1] - self.yz_values[i,1])
                            * (self.yz_values[i + 1, 0] - self.yz_values[i,0])
                            + self.yz_values[i,0],
                            h,
                        )
                    )
                    continue
            
            

        return np.array(yz_coords_new)
    
    def _calculate_wet_perimeter(self, h: float):
        perimeter = 0 
        for i in range(len(self.wet_profile)-1):
            if self.wet_profile[i,1] < h or self.wet_profile[i + 1,1] < h:
                perimeter += math.sqrt(( self.wet_profile[i,0] - self.wet_profile[i + 1, 0] ) ** 2 
                                         + ( self.wet_profile[i,1] - self.wet_profile[i + 1, 1] ) ** 2)
        return perimeter

    def _calc_froude(self):
        if self.wet_area > 0:
            return self.discharge/(self.wet_area  *math.sqrt(9.81 * max(self.water_level - min(self.yz_values[:, 1]),0)))
        else:
            return 0

    def _calc_freeboard_left(self):
        return self.yz_values[0, 1] - self.water_level

    def _calc_freeboard_right(self):
        return self.yz_values[-1, 1] - self.water_level     

    def calc_discharge(self, water_level):
        self._reset_calculation_output()
        self.water_level = water_level
        if water_level<=min(self.yz_values[:, 1]):
            self.discharge = 0
        else:
            if self.b_b:
                self.k_bb=self.bos_en_bijkerk_to_k_manning(self.k_manning, max(self.water_level - min(self.yz_values[:, 1]), 0))
                self.discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            else:
                self.discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            
        return self.discharge

    def calc_water_level(self, discharge):
        """
        Iteratively determines which water level occurs at a given discharge. First, the discharge at a water depth of
        1 meter is determined. The water level is then either increased or decreased until the right discharge is found.
        """
        self._reset_calculation_output()
        self.discharge=discharge
        self.water_level=min(self.yz_values[:, 1]) + 1
        delta=0.1
        if self.b_b:
            self.k_bb=self.bos_en_bijkerk_to_k_manning(self.k_manning, max(self.water_level - min(self.yz_values[:, 1]), 0))
            check_discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
        else:
            check_discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius)
        if check_discharge  >discharge:
            delta = -delta
            while check_discharge > discharge:
                self.water_level += delta
                if self.b_b:
                    self.k_bb = self.bos_en_bijkerk_to_k_manning(self.k_manning, max(self.water_level - min(self.yz_values[:, 1]), 0))
                    check_discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
                else:
                    check_discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            delta = -delta * 0.1
            while check_discharge < discharge:
                self.water_level += delta
                if self.b_b:
                    self.k_bb = self.bos_en_bijkerk_to_k_manning(self.k_manning, max(self.water_level - min(self.yz_values[:, 1]), 0))
                    check_discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
                else:
                    check_discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            delta = -delta  * 0.1
            while check_discharge  > discharge:
                self.water_level += delta
                if self.b_b:
                    self.k_bb = self.bos_en_bijkerk_to_k_manning(self.k_manning, max(self.water_level - min(self.yz_values[:, 1]), 0))
                    check_discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
                else:
                    check_discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            self.water_level -= delta * 0.5
            if self.b_b:
                self.k_bb = self.bos_en_bijkerk_to_k_manning(self.k_manning,max(self.water_level - min(self.yz_values[:,1]),0))
                check_discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            else:
                check_discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius) 
            if check_discharge > discharge:
                self.water_level += delta * 0.5
            else:
                self.water_level-=delta * 0.5
        elif check_discharge<discharge:
            while check_discharge < discharge:
                self.water_level += delta
                if self.b_b:
                    self.k_bb = self.bos_en_bijkerk_to_k_manning(self.k_manning,max(self.water_level - min(self.yz_values[:,1]),0))
                    check_discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
                else:
                    check_discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            delta = -delta * 0.1
            while check_discharge > discharge:
                self.water_level += delta
                if self.b_b:
                    self.k_bb=self.bos_en_bijkerk_to_k_manning(self.k_manning, max(self.water_level - min(self.yz_values[:, 1]), 0))
                    check_discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
                else:
                    check_discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            delta = -delta * 0.1
            while check_discharge < discharge:
                self.water_level += delta
                if self.b_b:
                    self.k_bb=self.bos_en_bijkerk_to_k_manning(self.k_manning, max(self.water_level - min(self.yz_values[:, 1]), 0))
                    check_discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
                else:
                    check_discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            self.water_level -= delta * 0.5
            if self.b_b:
                self.k_bb = self.bos_en_bijkerk_to_k_manning(self.k_manning, max(self.water_level - min(self.yz_values[:, 1]), 0))
                check_discharge = self._calc_discharge(self.k_bb, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            else:
                check_discharge = self._calc_discharge(self.k_manning, self.bottom_slope, self.wet_area, self.hydraulic_radius)
            if check_discharge < discharge:
                self.water_level += delta*0.5
            else:
                self.water_level -= delta*0.5                    
        return self.water_level 

    def _reset_calculation_output(self):
        self.water_level = None
        self.discharge = None    

if __name__ == "__main__":

    pf = Trapezoidal(30, 0.0008, 1, 1, 1, 0, 1)
    q = pf.calc_discharge(2)
    print(q)

    cs = Trapezoidal(
        k_manning=30,
        bottom_slope=0.0008,
        slope_left=1,
        slope_right=1,
        bed_width=1,
        bed_level=0,
        surface_level=1)
    cs.calc_discharge(1)
