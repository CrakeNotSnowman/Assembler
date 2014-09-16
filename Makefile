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
	rm -rf build
	python setup2.py build
	sudo python setup2.py install
	python setup3cpp.py build
	sudo python setup3cpp.py install

clean:
	rm -rf build
