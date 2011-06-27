### pvt.hcobj.R  (2011-06-22)
###
###    Compute empirical HC score from p-values
###
### Copyright 2009-11 Korbinian Strimmer
###
###
### This file is part of the `sda' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA



# private function (not documented)

# empirical HC objective from p-values
pvt.hcobj = function(pval)
{
  d = length(pval)
  #F = ecdf(pval)(pval)
  F = rank(pval, ties.method="max")/d
  
  v = F*(1-F)/d # variance
  v[v==0] = min( v[v > 0] ) # just to make sure we have no zero variance
  HC = abs(F-pval)/sqrt(v)

  return( HC )
}


