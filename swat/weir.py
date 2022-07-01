import numpy as np

class Weir():
    def __init__(self, crest: float, width: float, Cs: float, Cw: float, Cl: float, 
                 water_level_up: float, water_level_down: float, gate=np.inf, Cc=0.63):
        """
        Calculate discharge through a weir or orofice structure. determines free or submerged flow automatically
        crest = crest level of the structure (m above datum)
        width = width of the structure (m)
        Cs = loss coeffcietn of structure (this is used during submerged flow)
        Cw = weir coeffcient depends on the shape of the crest (e.g. sharp, broad, etc)
        Cl = lateral contraction coeffcient will constrict the flow width
        water_level_up = upstream water level (m above datum)
        water_level_down = downstream water level (m above datum)
        gate = gate height above crest (m)
        Cc = contraction coefficient, will constrain the flow height under the gate
        """

        self.crest = crest
        self.width = width
        self.gate = gate
        self.Cc = Cc
        self.Cw = Cw
        self.Cl = Cl
        self.Cs = Cs

        #check which water level is higher (determine positive of negative flow)
        if water_level_up > water_level_down:
            self.water_level1 = water_level_up
            self.water_level2 = water_level_down
            self.sign = 1
        else:
            self.water_level1 = water_level_down
            self.water_level2 = water_level_up
            self.sign = -1 
 
    @property
    def flow_check(self):        
        """
        Free or submerged flow
        """
        return self._flow_check()

    @property
    def velocity(self):
        """
        Flow velocity over crest
        """
        return self._velocity()
        
    @property
    def wet_area(self):
        """
        Wet area over crest
        """
        return self._wet_area()

    @property
    def discharge(self):
        """
        Calculated discharge
        """
        return self._discharge()

    @property
    def height_above_crest(self):
        """
        height above crest or opening height if gate is relevant
        """
        return self._height_above_crest()        

    def _flow_check(self):
        #check flow situation
        if self.water_level1 - self.crest < 1.5 * self.gate:
            # pure weir flow (gate is not relevant)
            if self.water_level1 - self.crest > 1.5 * ( self.water_level2 - self.crest) :
                #free flow
                flow_check = True
            else:
                #submerged flow
                flow_check = False
        else:
            # orifice flow (gate matters)
            if self.water_level2 <= self.crest+self.gate:
                #free flow
                flow_check = True
            else:
                #submerged flow
                flow_check = False                     
        return flow_check

    def _velocity(self):
        if self.water_level1 - self.crest < 1.5 * self.gate:
            if self.water_level1 > self.crest:
                if self.flow_check:
                    velocity = self.Cw*np.sqrt(2/3 * 9.81 * (self.water_level1 - self.crest)) * self.sign
                else:
                    velocity = self.Cs*np.sqrt(2 * 9.81 * (self.water_level1 - self.water_level2)) * self.sign
            else:
                velocity=0
        else:
            if self.gate>0:
                if self.flow_check:
                    velocity  =np.sqrt(2 * 9.81 * ( self.water_level1 - ( self.crest + self.Cc * self.gate ))) * self.sign
                else:
                    velocity = self.Cs*np.sqrt(2 * 9.81 * ( self.water_level1 - self.water_level2)) * self.sign
            else:
                velocity = 0
        return velocity

    def _wet_area(self):
        if self.water_level1-self.crest < 1.5 * self.gate:
            if self.flow_check:
                wet_area=max(self.Cl  *self.width * 2/3  *( self.water_level1 - self.crest ), 0)
            else:
                wet_area=max(self.Cl * self.width * min( self.water_level1 - self.crest - (self.velocity ** 2) / (2 * 9.81), self.gate) ,0)
        else:
            wet_area=max(self.Cl * self.width * self.Cc * self.gate, 0)
        
        return wet_area
        
    def _discharge(self): 
        discharge=self.wet_area * self.velocity
        return discharge    

    def _height_above_crest(self):
        if self.water_level1 - self.crest < 1.5 * self.gate:
            return max(self.water_level1 - self.crest, 0)
        else:
            return self.gate    

if __name__ == "__main__":
    struc = Weir(
                crest=2,
                width=3,
                gate=0.5,
                Cc=0.63,
                Cw=1,
                Cs=0.9,
                Cl=1,
                water_level_up=3.5,
                water_level_down=2)

    print(struc.flow_check) 
    print(struc.discharge) 
    print(struc.wet_area) 
    print(struc.velocity) 
