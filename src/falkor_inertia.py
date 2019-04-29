from sympy import symbols
from sympy.physics.mechanics import ReferenceFrame, Point, inertia, RigidBody


class RigidBody2(RigidBody):

    def parallel_axis(self, point):
        """Returns the inertia of the body about another point."""
        # TODO : What if the new point is not fixed in the rigid body's frame?
        a, b, c = self.masscenter.pos_from(point).to_matrix(self.frame)
        return self.central_inertia + self.mass * inertia(self.frame,
                                                          b**2 + c**2,
                                                          c**2 + a**2,
                                                          a**2 + b**2,
                                                          -a * b,
                                                          -b * c,
                                                          -a * c)


# this section is trying to do the inertia addition with a composite head
#h_u, h_l, l_u, l_l, w = symbols('h_u, h_l, l_u, l_l, w')
#x_t, y_t, z_t = symbols('x_t, y_t, z_t')
#m_u, m_l, m_t = symbols('m_u, m_l, m_t')
#
#F = ReferenceFrame('F')
#
#I_u = m_u / 12 * inertia(F, h_u**2 + w**2, h_u**2 + l_u**2, l_u**2 + w**2)
#I_l = m_l / 12 * inertia(F, h_l**2 + w**2, h_l**2 + l_l**2, l_l**2 + w**2)
#
#o = Point('o')
#uo = o.locatenew('u_o', l_u / 2 * F.x + w / 2 * F.y + h_u / 2 * F.z)
#lo = o.locatenew('l_o', l_l / 2 * F.x + w / 2 * F.y + h_l / 2 * F.z)
#to = o.locatenew('t_o', (m_u * uo.pos_from(o) + m_l * lo.pos_from(o)) / m_t)
#
#U = RigidBody2('U', uo, F, m_u, (I_u, uo))
#L = RigidBody2('U', lo, F, m_l, (I_l, lo))
#
#I_t = U.parallel_axis(to) + L.parallel_axis(to)

# F: earth reference frame
# H: front assembly (handlebar and fork)
# D: dragon head
# T: combined dragon head and fron assembly

# this models the head simply as a single cuboid
inch_per_meter = 0.0254
ipm = inch_per_meter

w = symbols('w')
lz, lx, ly = symbols('lz, lx, ly')  # dimensions of falkor head cuboid
xH, yH, zH = symbols('x_H, y_H, z_H')  # front assembly CoM location
xT, yT, zT = symbols('x_T, y_T, z_T')  # dragon + front assembly CoM location
mD, mH, mT = symbols('m_D, m_H, m_T')
IHxx, IHyy, IHzz, IHxz = symbols('I_Hxx, I_Hyy, I_Hzz, I_Hxz')

F = ReferenceFrame('F')

IH = inertia(F, IHxx, IHyy, IHzz, izx=IHxz)  # central inertia of front assembly
ID = mD / 12 * inertia(F, lz**2 + ly**2, lz**2 + lx**2, lx**2 + ly**2)  # central inertia of dragon head

o = Point('o')  # bicycle origin (contact of rear wheel)
# locate the mass center of the handlebar assembly with respect to origin
Ho = o.locatenew('Ho', xH * F.x + yH * F.y + zH * F.z)
# locate the mass center of the dragon head  with respect to origin
Do = o.locatenew('D_o', (w + 7.5 * ipm) * F.x - (27.5 + 9.0) * ipm * F.z)
To = o.locatenew('T_o', (mD * Do.pos_from(o) + mH * Ho.pos_from(o)) / (mD + mH))

H = RigidBody2('H', Ho, F, mH, (IH, Ho))
D = RigidBody2('D', Do, F, mD, (ID, Do))

IT = H.parallel_axis(To) + D.parallel_axis(To)

par = {w: 1.02,
       ly: 22 * ipm,
       lx: 36 * ipm,
       lz: 9 * ipm,
       mD: 5.987419,
       mH: 4,
       xH: 0.9,
       yH: 0.0,
       zH: -0.7,
       IHxx: 0.05892,
       IHyy: 0.06,
       IHzz: 0.00708,
       IHxz: -0.00756}
