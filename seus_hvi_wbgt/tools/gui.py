#!/usr/bin/env python3

import logging

import numpy

from metpy.units import units

from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QFileDialog, QLabel, QPushButton, QComboBox, QTextEdit, QLineEdit
from PyQt5.QtWidgets import QVBoxLayout, QHBoxLayout, QGridLayout, QAction
from PyQt5.QtGui import QIcon, QPixmap
from PyQt5.QtCore import Qt

from seus_hvi_wbgt.wbgt import wbgt

class WBGT_Gui( QMainWindow ):

  def __init__(self, *args, **kwargs):

    super().__init__(*args, **kwargs)

    mainWidget = QWidget()
    layout = QGridLayout()

    monthLabel     = QLabel( 'Month' )
    dayLabel       = QLabel( 'Day' )
    latLabel       = QLabel( 'Latitude' )
    solarLabel     = QLabel( 'Solar (W/m^2)' )
    presLabel      = QLabel( 'Pressure' )
    fcstTempLabel  = QLabel( 'Fcst Max Temp (F)' )
    dewTempLabel   = QLabel( 'DewPoint (F)' )
    windLabel      = QLabel( 'Wind Speed (mph)' )
    cloudLabel     = QLabel( 'Cloud Cover (%)' )

    liljegrenLabel = QLabel( 'Liljegren' )
    dimiceliLabel  = QLabel( 'Dimiceli' )
    bernardLabel   = QLabel( 'Bernard' )

    self.month     = QLineEdit()
    self.day       = QLineEdit()
    self.lat       = QLineEdit()
    self.solar     = QLineEdit()
    self.pres      = QLineEdit()
    self.fcstTemp  = QLineEdit()
    self.dewTemp   = QLineEdit()
    self.wind      = QLineEdit()
    self.cloud     = QLineEdit()

    self.presUnits = QComboBox()
    self.presUnits.addItems( ['hPa', 'inHg'] )
    self.presUnits.currentIndexChanged.connect( self.valChanged )

    self.month.setText(  '1')
    self.day.setText(    '1')
    self.lat.setText(    '0')
    self.solar.setText(  '1000')
    self.pres.setText(   '1000')
    self.fcstTemp.setText( '80' )
    self.dewTemp.setText(  '60' )
    self.wind.setText(      '5' )
    self.cloud.setText(     '0' )

    liljegren      = QLabel( '' )
    dimiceli       = QLabel( '' )
    bernard        = QLabel( '' )

    self.month.textChanged.connect(    self.valChanged )
    self.day.textChanged.connect(      self.valChanged ) 
    self.lat.textChanged.connect(      self.valChanged ) 
    self.solar.textChanged.connect(    self.valChanged ) 
    self.pres.textChanged.connect(     self.valChanged ) 
    self.fcstTemp.textChanged.connect( self.valChanged ) 
    self.dewTemp.textChanged.connect(  self.valChanged ) 
    self.wind.textChanged.connect(     self.valChanged ) 
    self.cloud.textChanged.connect(    self.valChanged ) 

    layout.addWidget( monthLabel,      0, 0) 
    layout.addWidget( dayLabel,        1, 0) 
    layout.addWidget( latLabel,        2, 0) 
    layout.addWidget( solarLabel,      3, 0) 
    layout.addWidget( presLabel,       4, 0) 
    layout.addWidget( fcstTempLabel,   5, 0) 
    layout.addWidget( dewTempLabel,    6, 0) 
    layout.addWidget( windLabel,       7, 0) 
    layout.addWidget( cloudLabel,      8, 0) 
    layout.addWidget( liljegrenLabel,  9, 0) 
    layout.addWidget( dimiceliLabel,  10, 0) 
    layout.addWidget( bernardLabel,   11, 0) 

    layout.addWidget( self.month,      0, 1) 
    layout.addWidget( self.day,        1, 1) 
    layout.addWidget( self.lat,        2, 1) 
    layout.addWidget( self.solar,      3, 1) 
    layout.addWidget( self.pres,       4, 1) 
    layout.addWidget( self.fcstTemp,   5, 1) 
    layout.addWidget( self.dewTemp,    6, 1) 
    layout.addWidget( self.wind,       7, 1) 
    layout.addWidget( self.cloud,      8, 1) 
    layout.addWidget( liljegren,       9, 1) 
    layout.addWidget( dimiceli,       10, 1) 
    layout.addWidget( bernard,        11, 1) 

    layout.addWidget( self.presUnits,  4, 2)

    self.wbgt = {
      'liljegren' : liljegren,
      'dimiceli'  : dimiceli,
      'bernard'   : bernard   
    }

    mainWidget.setLayout( layout )
    self.setCentralWidget( mainWidget )

    self.valChanged()

    self.show()


  def valChanged( self, *args, **kwargs ):

    try:
      month    = numpy.asarray( [  int( self.month.text(     ) )] )
      day      = numpy.asarray( [  int( self.day.text(       ) )] )
      lat      = numpy.asarray( [float( self.lat.text(       ) )] ) 
      solar    = numpy.asarray( [float( self.solar.text(     ) )] ) * units('watt/m**2')
      pres     = numpy.asarray( [float( self.pres.text(      ) )] )
      fcstTemp = numpy.asarray( [float( self.fcstTemp.text(  ) )] ) * units.degree_Fahrenheit
      dewTemp  = numpy.asarray( [float( self.dewTemp.text(   ) )] ) * units.degree_Fahrenheit
      wind     = numpy.asarray( [float( self.wind.text(      ) )] ) * units.mph
      cloud    = numpy.asarray( [float( self.cloud.text(     ) )] ) 
    except:
      return

    year   = numpy.asarray( [2000])
    hour   = numpy.asarray( [  12])
    minute = numpy.asarray( [   0])
    lon    = numpy.asarray( [   0.0])
    pres   = pres * units( self.presUnits.currentText() )
    print( pres )

    for method, widget in self.wbgt.items():
      Twbg = wbgt( method, 
        lat, lat, year, month, day, hour, minute, solar, pres, fcstTemp, dewTemp, wind
      )
      widget.setText( str( Twbg['Twbg'][0] * 9.0/5.0 + 32 ) )


