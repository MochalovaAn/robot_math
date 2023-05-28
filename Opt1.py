import Parts.Constans as co
import Parts.Visual as vi
import Parts.MainFunction as mf
import Parts.differ as di


print(co.L)
print(((co.m+co.M)*co.L/(co.m*co.R))**2)
# Sist = [co.psi0, co.w0, co.ksi10 ,co.ksi20]
#
# def eevent(t,Sist): return mf.getY(t, Sist)/(mf.getX(t, Sist[1]) + 0.0000001) + 10000000000
# step = 0.001
#
# res = di.resolveSist(mf.SoprSist, 50, Sist, eevent, step)
#
# vi.plotPsi(res)