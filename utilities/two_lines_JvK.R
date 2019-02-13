# These scripts reproduce the analysis in the paper: van Kempen et al.,
# (2018) 'Behavioural and neural signatures of perceptual evidence
# accumulation are modulated by pupil-linked arousal'
# 
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the 
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#   
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# 
# Jochem van Kempen, 2018
# Jochemvankempen@gmail.com
# https://github.com/jochemvankempen/2018_Monash
# 
# -------------------------------------------------------------------------
# Code to test whether there is a U-shaped relationship.
# This code was inspired by the code written by Uri Simonsohn for which
# the exact code can be found at http://webstimate.org/twolines/
#
# It is adapted to be able to apply lmer models in order to be able to use random intercepts etc.
# Here I apply it to averaged (binned) data, which simplifies the procedure substantially 

reg2_JvK=function(x,y,xc,data,excludeXc)
{
  #Create new variables
  if (excludeXc==0)
    {xlow1=ifelse(x<=xc,x-xc,0)}     #xlow=x-xc when x<xc, 0 otherwise
  else 
    {xlow1=ifelse(x<xc,x-xc,0)}
  xhigh1=ifelse(x>xc,x-xc,0)     #xhigh=x when x<xmax, 0 otherwise
  high1=ifelse(x>xc,1,0)         #high dummy, allows interruption
  
  #Run the regressions 
  lm1=lmerTest::lmer(y~xlow1+xhigh1+high1+(1|Subject),
                 data = data, na.action = na.omit, REML=FALSE)     #estimate regression
  lmc1=summary(lm1)$coefficients   #Get b,t,se,

  #Turn into individual results (scalars)
  b1=lmc1[2,1]
  p1=lmc1[2,5]
  
  b2=lmc1[3,1]
  p2=lmc1[3,5]
  
  #Is the u-shape significant?
  usig =ifelse(b1*b2<0 & p1<.0167 & p2<.0167,1,0)                     
  #All results
  res = usig
  res
  
}