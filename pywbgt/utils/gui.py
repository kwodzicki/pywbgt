"""
GUI for computing WBGT

"""

import numpy

from metpy.units import units

from PyQt5.QtWidgets import (
    QMainWindow,
    QWidget,
    QLabel,
    QComboBox,
    QLineEdit,
    QGridLayout,
)

from seus_hvi_wbgt.wbgt import wbgt

class WBGTGui( QMainWindow ):
    """
    Create GUI for WBGT

    """

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        main_widget = QWidget()
        layout = QGridLayout()

        month_label     = QLabel( 'Month' )
        day_label       = QLabel( 'Day' )
        lat_label       = QLabel( 'Latitude' )
        solar_label     = QLabel( 'Solar (W/m^2)' )
        pres_label      = QLabel( 'Pressure' )
        fcst_temp_label = QLabel( 'Fcst Max Temp (F)' )
        dew_temp_label  = QLabel( 'DewPoint (F)' )
        wind_label      = QLabel( 'Wind Speed (mph)' )
        cloud_label     = QLabel( 'Cloud Cover (%)' )

        liljegren_label = QLabel( 'Liljegren' )
        dimiceli_label  = QLabel( 'Dimiceli' )
        bernard_label   = QLabel( 'Bernard' )

        self.month     = QLineEdit()
        self.day       = QLineEdit()
        self.lat       = QLineEdit()
        self.solar     = QLineEdit()
        self.pres      = QLineEdit()
        self.fcst_temp = QLineEdit()
        self.dew_temp  = QLineEdit()
        self.wind      = QLineEdit()
        self.cloud     = QLineEdit()

        self.pres_units = QComboBox()
        self.pres_units.addItems( ['hPa', 'inHg'] )
        self.pres_units.currentIndexChanged.connect( self.val_changed )

        self.month.setText(  '1')
        self.day.setText(    '1')
        self.lat.setText(    '0')
        self.solar.setText(  '1000')
        self.pres.setText(   '1000')
        self.fcst_temp.setText( '80' )
        self.dew_temp.setText(  '60' )
        self.wind.setText(      '5' )
        self.cloud.setText(     '0' )

        liljegren      = QLabel( '' )
        dimiceli       = QLabel( '' )
        bernard        = QLabel( '' )

        self.month.textChanged.connect(     self.val_changed )
        self.day.textChanged.connect(       self.val_changed )
        self.lat.textChanged.connect(       self.val_changed )
        self.solar.textChanged.connect(     self.val_changed )
        self.pres.textChanged.connect(      self.val_changed )
        self.fcst_temp.textChanged.connect( self.val_changed )
        self.dew_temp.textChanged.connect(  self.val_changed )
        self.wind.textChanged.connect(      self.val_changed )
        self.cloud.textChanged.connect(     self.val_changed )

        layout.addWidget( month_label,      0, 0)
        layout.addWidget( day_label,        1, 0)
        layout.addWidget( lat_label,        2, 0)
        layout.addWidget( solar_label,      3, 0)
        layout.addWidget( pres_label,       4, 0)
        layout.addWidget( fcst_temp_label,   5, 0)
        layout.addWidget( dew_temp_label,    6, 0)
        layout.addWidget( wind_label,       7, 0)
        layout.addWidget( cloud_label,      8, 0)
        layout.addWidget( liljegren_label,  9, 0)
        layout.addWidget( dimiceli_label,  10, 0)
        layout.addWidget( bernard_label,   11, 0)

        layout.addWidget( self.month,      0, 1)
        layout.addWidget( self.day,        1, 1)
        layout.addWidget( self.lat,        2, 1)
        layout.addWidget( self.solar,      3, 1)
        layout.addWidget( self.pres,       4, 1)
        layout.addWidget( self.fcst_temp,   5, 1)
        layout.addWidget( self.dew_temp,    6, 1)
        layout.addWidget( self.wind,       7, 1)
        layout.addWidget( self.cloud,      8, 1)
        layout.addWidget( liljegren,       9, 1)
        layout.addWidget( dimiceli,       10, 1)
        layout.addWidget( bernard,        11, 1)

        layout.addWidget( self.pres_units,  4, 2)

        self.wbgt = {
          'liljegren' : liljegren,
          'dimiceli'  : dimiceli,
          'bernard'   : bernard,
        }

        main_widget.setLayout( layout )
        self.setCentralWidget( main_widget )

        self.val_changed()

        self.show()


    def val_changed( self, *args, **kwargs ):
        """
        Run when value in text box changes

        """

        try:
            month     = numpy.asarray( [  int( self.month.text(     ) )] )
            day       = numpy.asarray( [  int( self.day.text(       ) )] )
            lat       = numpy.asarray( [float( self.lat.text(       ) )] )
            solar     = numpy.asarray( [float( self.solar.text(     ) )] ) * units('watt/m**2')
            pres      = numpy.asarray( [float( self.pres.text(      ) )] )
            fcst_temp = numpy.asarray( [float( self.fcst_temp.text( ) )] ) * units.degree_Fahrenheit
            dew_temp  = numpy.asarray( [float( self.dew_temp.text(  ) )] ) * units.degree_Fahrenheit
            wind      = numpy.asarray( [float( self.wind.text(      ) )] ) * units.mph
            #cloud     = numpy.asarray( [float( self.cloud.text(     ) )] )
        except:
            return

        year   = numpy.asarray( [2000])
        hour   = numpy.asarray( [  12])
        minute = numpy.asarray( [   0])
        lon    = numpy.asarray( [   0.0])
        pres   = pres * units( self.pres_units.currentText() )
        print( pres )

        for method, widget in self.wbgt.items():
            temp_wbg = wbgt( method,
                lat, lat, year, month, day, hour, minute, solar, pres, fcst_temp, dew_temp, wind
            )
            widget.setText( str( temp_wbg['Twbg'][0] * 9.0/5.0 + 32 ) )
