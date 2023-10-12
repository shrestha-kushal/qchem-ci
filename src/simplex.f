      subroutine simplex (fn, n, nn, start, xmin, ylo, step, ical, iter)
c 
c     Dr. John J. O'Dea            Whitten group
c 
c     November 24, 2003
c 
c     optimization using simplex formation
c     with Nelder-Mead rules
c 
c     References:
c 
c     Nelder, J. A., and Mead, R.,
c     "A Simplex Method for Function Minimization",
c     the computer journal, vol. 7, pages 308-313, (1965).
c 
c     Olsson, D. M.,
c     "A Sequential Simplex Program for Solving
c     Minimization Problems",
c     journal of quality technology, vol. 6, pages 53-57,
c     (1974).
c 
c      Largarias, J. C.; Reeds, J. A.;, Wright, M. H; et al.,
c      "Convergence Properties of the Nelder-Mead Simplex
c      Method in Low Dimensions", SIAM Journal on Optimization,
c      Vol. 9, pages 112-147, (1998).
c 
c 
c     'n' = number of dimensions in simplex space = number
c     of arguments in function 'fn'.  'nn' = 'n'+1 = number
c     of simplex vertices.
c 
c     'start' contains the coordinates of the initial
c     vertex, ie., the initial estimates of the variables
c     of 'fn'.
c 
c     'xmin' contains the optimum vertex coordinates
c     from the last iteration.
c 
c     'ylo' is the function value of the most favorable
c     vertex.
c 
c     'reqmin' is a convergence criterion.
c 
c     'step' is an array 'n' elements long which contains
c     the initial increments for the coordinates.
c 
c     'ical' counts the number of objective function calls
c     required to find the optimum.
c 
c     'itmax' is the maximum number of iterations allowed.
c 
c     'iter' counts the number of iterations performed.
c     entry into this subroutine is counted as one iteration
c     even if the first iteration is not completed.
c 
c     'x(n,nn)' is a two dimensional array which contains
c     'n' number of coordinates of 'nn' number of vertices,
c     ie., 'nn' designates the vertex which is specified by
c     'n' number of coordinates.
c 
c     'y' is an array 'nn' elements long which contains the
c     function values of 'nn' vertices.
c 
c     'xbar' is an array 'n' elements long which
c     contains the coordinates of the centroid of
c     the most favorable 'n' vertices.
c 
c     'xr' is an array 'n' elements long which contains
c     the coordinates of a vertex after reflection through
c     the centroid.
c 
c     'yr' is the function value at the
c     coordinates specified by 'xr'.
c 
c     'xnew' is an array 'n' elements long which
c     contains the coordinates of a new vertex
c     calculated by expansion or contraction.
c 
c     'ynew' is the function value at the
c     coordinates specified by 'xnew'.
c 
c     'kc' is the maximum number of function calls
c     allowed.
c 
      implicit none
      integer n,nn
      double precision start(n),step(n),xmin(n),
     *x(n,nn),xr(n),xnew(n),xbar(n),y(nn)
      double precision fn,alpha,beta,gama,dabit,reqmin,
     *temp,ylo,yhi,yoldlo,coord1,coord2,dchk,z,yr,ynew
      integer i,j,itmax,ical,kc,iter,ilo,ihi,l,ibest
      logical lan
c 
      external fn
c 
c     define contraction and expansion coefficients.
c 
      data alpha,beta,gama / 1.d0,0.5d0,2.d0 /
      data dabit / 2.0407d-35 /
      data itmax / 40000 /
      data kc / 90000 /
      data reqmin / 1.0d-7 /
c 
      ical=0
      iter=1
c 
c     'nn' is the number of vertices in the simplex.
c 
c      nn=n+1
c 
c     construct the initial simplex.
c     use the corner algorithm, ie., increment the
c     coordinates of the initial estimate space in
c     one direction at a time by the specfied step.
c 
c     define vertex 'nn' using the coordinates
c     of the initial estimates.
c 
      do 110 i=1,n
      x(i,nn)=start(i)
110   continue
c 
c     find the objective function for the initial estimate.
c 
      y(nn)=fn(start)
      ical=ical+1
c 
c     the coordinates of the remaining vertices
c     are defined by adding a step to each respective
c     dimension.
c 
      do 130 j=1,n
      temp=start(j)
      start(j)=temp+step(j)
c 
      do 120 i=1,n
      x(i,j)=start(i)
120   continue
c 
c     find the objective function for the remaining vertices.
c 
      y(j)=fn(start)
      ical=ical+1
      start(j)=temp
130   continue
c 
c     the initial simplex construction is complete.
c     begin simplex formation using Nelder-Mead rules.
c 
c     find highest and lowest objective function values.
c 
140   ylo=y(1)
      yhi=y(1)
      ilo=1
      ihi=1
c 
      do 160 i=2,nn
      if (y(i).ge.ylo) go to 150
      ylo=y(i)
      ilo=i
150   if (y(i).le.yhi) go to 160
      yhi=y(i)
      ihi=i
160   continue
c 
      do 800 i=1,nn
 800  write(18,*)y(i)
      call flush(18)
c     ylo = y(ilo) is the smallest function value of
c     all the vertices.
c 
c     yhi = y(ihi) is the largest function value of
c     all the vertices.
c 
      if (ical.le.nn) yoldlo=ylo
      if (ical.le.nn) go to 170
      if (ylo.ge.yoldlo) go to 170
      yoldlo=ylo
c 
c     stop if the maximum number of iterations is reached.
c 
      iter=iter+1
      if (iter.ge.itmax) go to 370
c 
c     check for existence of 'optim.run'.
c     if this file does not exist, stop iterations,
c     and return results.
c 
c170   inquire (file='optim.run',exist=lan)
c     if (.not.lan) go to 370
 170  continue
c 
c     check for convergence of coordinates.
c 
      do 200 i=1,n
      coord1=x(i,1)
      coord2=coord1
c 
      do 190 j=2,nn
      if (x(i,j).ge.coord1) go to 180
      coord1=x(i,j)
180   if (x(i,j).le.coord2) go to 190
      coord2=x(i,j)
190   continue
c 
      dchk=(coord2+dabit)/(coord1+dabit)-1.d0
      if (dabs(dchk).gt.reqmin) go to 210
200   continue
c 
c     quit, convergence has been achieved.
c 
      go to 370
c 
c     quit if the maximum number of function calls is reached.
c 
210   if (ical.ge.kc) go to 370
c 
c     now calculate 'xbar', the centroid of all the simplex
c     vertices except one. the vertex not included has
c     the largest objective function value.
c 
      do 230 i=1,n
      z=0.d0
c 
      do 220 j=1,nn
      z=z+x(i,j)
220   continue
c 
      z=z-x(i,ihi)
      xbar(i)=z/n
230   continue
c 
c     find the coordinates of the reflection of x(i,ihi)
c     through the centroid, ie., find 'xr'.
c 
      do 240 i=1,n
      xr(i)=(1.d0+alpha)*xbar(i)-alpha*x(i,ihi)
240   continue
c 
c     find the function at the reflected coordinates, 'xr'.
c 
      yr=fn(xr)
      ical=ical+1
c 
c     test the new function value, 'yr', against the
c     most favorable vertex, 'ylo', and invoke the rules
c     for new simplex formation.
c 
      if (yr.ge.ylo) go to 280
      if (ical.ge.kc) go to 350
c 
c     the reflection produced a new minimum,
c     so perform an extension.
c 
      do 250 i=1,n
      xnew(i)=gama*xr(i)+(1.0-gama)*xbar(i)
250   continue
c 
      ynew=fn(xnew)
      ical=ical+1
      if (ynew.ge.yr) go to 350
c 
c     accept the extension or contraction, ie., replace
c     the vertex x(i,ihi) with 'xnew'.
c 
260   do 270 i=1,n
      x(i,ihi)=xnew(i)
270   continue
      y(ihi)=ynew
c 
c     start another iteration.
c 
      go to 140
c 
c     no extension is indicated. test 'yr' against the
c     other vertices.
c 
280   l=0
c 
      do 290 i=1,nn
      if (y(i).gt.yr) l=l+1
290   continue
c 
c     if more than one vertex is greater than 'yr'
c     accept the reflected vertex.
c 
      if (l.gt.1) go to 350
c 
c     if no vertex is greater than 'yr', contraction
c     occurs on the y(ihi) side of the centroid.
c 
      if (l.eq.0) go to 310
c 
c     if only one vertex is greater than 'yr', contraction
c     occurs on the reflection side of the centroid.
c 
      do 300 i=1,n
      x(i,ihi)=xr(i)
300   continue
c 
      y(ihi)=yr
c 
c     contraction on the y(ihi) side of the centroid.
c 
310   if (ical.ge.kc) go to 370
c 
      do 320 i=1,n
      xnew(i)=beta*x(i,ihi)+(1.0-beta)*xbar(i)
320   continue
      ynew=fn(xnew)
      ical=ical+1
      if (ynew.lt.y(ihi)) go to 260
c 
c     contract the whole simplex.
c 
      do 340 j=1,nn
      do 330 i=1,n
      x(i,j)=(x(i,j)+x(i,ilo))*0.5d0
      xmin(i)=x(i,j)
330   continue
      y(j)=fn(xmin)
340   continue
c 
      ical=ical+nn
      if (ical.ge.kc) go to 370
c 
c     start another iteration.
c 
      go to 140
c 
c     accept the reflected vertex.
c 
350   do 360 i=1,n
      x(i,ihi)=xr(i)
360   continue
      y(ihi)=yr
c 
c     go start another iteration.
c 
      go to 140
c 
c     the optimization is finished so
c     save the optimum coordinates.
c 
370   do 380 i=1,n
      xmin(i)=x(i,nn)
380   continue
c 
c     compare 'nn' objective function values and
c     identify the best (most favorable) vertex.
c 
      ylo=y(1)
      ibest=1
c 
      do 390 j=2,nn
      if (y(j).ge.ylo) go to 390
      ylo=y(j)
      ibest=j
390   continue
c 
c     save the coordinates of the best vertex.
c 
      do 400 i=1,n
      xmin(i)=x(i,ibest)
400   continue
c 
c     finally, evaluate the objective function
c     at the optimal coordinates and return.
c 
      ylo=fn(xmin)
      ical=ical+1
c 
      return
      end


