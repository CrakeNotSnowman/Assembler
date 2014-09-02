#/*****************************************************************/
#/*  University of Nebraska-Lincoln                               */
#/*  Department of Electrical Engineering                         */
#/*  Adam "I don't know what I'm doing" Burbach                   */
#/*  			and				          */
#/*  Keith "Has joined the chaos" Murray			  */
#/*		     featuring					  */
#/*  Jacob "You owe me a doughnut" Bohac			  */	
#/*****************************************************************/


all:
	python setup2.py build
	sudo python setup2.py install

clean:
	rm -rf build
