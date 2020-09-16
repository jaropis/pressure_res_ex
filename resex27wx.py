from glob import glob
import os
import wx.lib.filebrowsebutton as filebrowse
import wx
import matplotlib

matplotlib.use("WXAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import (
    FigureCanvasWxAgg as FigCanvas,
    NavigationToolbar2WxAgg as NavigationToolbar,
)
from pylab import plot, show, legend, fill_between, subplot, xlabel, ylabel, savefig
from matplotlib.patches import FancyArrowPatch
from numpy import arange, array, diff, exp
from scipy import where
from scipy.optimize import leastsq, fmin
from scipy.integrate import simps, trapz
from csv import reader


class ListaDwustronna:
    """ta klasa jest mi potrzebna do guzikow previous next,
    ktore beda przechodzily przez liste plikow"""

    def __init__(self, lista):
        self.iterator = 0
        self.lista = lista
        self.dlugosc = len(self.lista)

    def get_first(self):
        return self.lista[0]

    def get_current(self):
        return self.lista[self.iterator]

    def rewind_to_position(self, filename):
        # w tej metodzie wyjdzie, ze jest zly typ pliku, wiec wyjatek tego typu (IndexError) zostanie zlapany
        print([i for i, x in enumerate(self.lista) if x == filename])
        try:
            self.iterator = [i for i, x in enumerate(self.lista) if x == filename][
                0
            ]  # w tym wyrazeniu listowym szukam, ktory numer ma dany plik
        except IndexError:
            app.frame.rysuj_on_error()
            return False
        return True

    def get_next(self):
        if self.iterator + 1 < self.dlugosc:
            self.iterator += 1
        else:
            self.iterator = 0
        return self.lista[self.iterator]

    def get_previous(self):
        if self.iterator - 1 >= 0:
            self.iterator -= 1
        else:
            self.iterator = self.dlugosc - 1
        return self.lista[self.iterator]

    def drukuj(self):
        """this method is only for testing"""
        print(self.lista)


class ReservoirFrame(wx.Frame):
    title = "ReserEx"

    def __init__(self):
        self.lista_plikow = ListaDwustronna([])
        self.dirname = ""
        self.filename = ""  # ta zmienna bedzie trzymala nazwe biezacego pliku (do uzycia w slowniku trzymajacym wyniki)
        self.diastolemethod = "Second Deriv"
        self.falloff = "full"
        self.czy_ruszane = False  # ta zmienna trzyma informacje o tym, czy jakikolwiek katalog byl otwarty. Jezeli nie, to nie bedzie zapisywany slownik wynikow. Do tej pory powodowalo to problemy, bo jezeli wylaczylo sie program bez zrobienia czegokolwiek, to zapisywany byl pusty slownik, ktory nadpisywal istniejacy slownik wynikowy i corruptowal caly program
        self.nodiast = False  ## czy uzyc algorytmu z prezentacji
        self.result_file = "nazwa"  # tu bede trzymal wyniki
        self.show_diastole = (
            False  ## czy rysowac punkt odpowiadajacy poczatkowi rozkurczu
        )
        self.dwietrzecie = False  ## czy uzywac 2/3 fall-offu czy calosci diastoli
        self.diastmethod = [
            True,
            False,
            False,
        ]  ## 1 to bedzie second derivative, 2 dicrotic notch, 3 sphygmocor
        self.results_dictionary = {}
        wx.Frame.__init__(self, None, -1, self.title)
        # self.create_status_bar()
        self.create_main_panel()
        self.create_menu()

    def get_pliki(self, sciezka):
        os.chdir(sciezka)
        return sorted(glob("*.txt"))

    def create_main_panel(self):
        self.panel = wx.Panel(self)
        global fig
        fig = plt.figure()
        # fig.subplots_adjust(left=0.03, bottom=0.03, top=0.97, right=0.98)
        self.canvas = FigCanvas(self.panel, -1, fig)
        # self.fig = Figure((5.0, 4.0), dpi=self.dpi)
        self.drawbuttonNext = wx.Button(self.panel, -1, "Next")
        self.drawbuttonNext.Bind(wx.EVT_BUTTON, self.on_Next)
        self.drawbuttonPrevious = wx.Button(self.panel, -1, "Previous")
        self.drawbuttonPrevious.Bind(
            wx.EVT_BUTTON, self.on_Previous, self.drawbuttonPrevious
        )

        self.drawbuttonRepeat = wx.Button(self.panel, -1, "Repeat")
        self.drawbuttonRepeat.Bind(wx.EVT_BUTTON, self.on_Repeat)

        self.drawbuttonCheckbox = wx.CheckBox(self.panel, -1, "Show diastole", (10, 10))
        self.drawbuttonCheckbox.Bind(wx.EVT_CHECKBOX, self.on_CheckBox)
        self.drawbuttonCheckbox2 = wx.CheckBox(
            self.panel, -1, "Use 2/3 of fall-off", (10, 10)
        )

        self.drawbuttonCheckboxNoD = wx.CheckBox(self.panel, -1, "No diast", (10, 10))

        self.drawradio = wx.RadioButton(
            self.panel, -1, label="Second Derivative", style=wx.RB_GROUP
        )
        self.drawradio2 = wx.RadioButton(self.panel, -1, label="D Notch")
        self.drawradio3 = wx.RadioButton(self.panel, -1, label="Sphygmocor")

        self.drawbuttonBatch = wx.Button(self.panel, -1, "Batch")
        self.drawbuttonBatch.Bind(wx.EVT_BUTTON, self.on_Batch)

        self.drawbuttonSave = wx.Button(self.panel, -1, "Save")
        self.drawbuttonSave.Bind(wx.EVT_BUTTON, self.on_Save)

        self.Bind(wx.EVT_CLOSE, self.on_Close)
        # self.toolbar = NavigationToolbar(self.canvas)## disabluje Toolbara na dole rysunku
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        # self.vbox.Add(self.toolbar, 0, wx.EXPAND)## disabluje Toolbara na dole rysunku
        self.vbox.AddSpacer(10)

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        # self.hbox.Add(self.textbox, 0, border=3, flag=flags)
        self.hbox.Add(self.drawbuttonPrevious, 0, border=3, flag=flags)
        self.hbox.Add(self.drawbuttonNext, 0, border=3, flag=flags)
        self.hbox.Add(self.drawbuttonRepeat, 0, border=3, flag=flags)
        self.hbox.Add(self.drawbuttonCheckbox, 0, border=3, flag=flags)
        self.hbox.Add(self.drawbuttonCheckbox2, 0, border=3, flag=flags)
        self.hbox.Add(self.drawbuttonCheckboxNoD, 0, border=3, flag=flags)
        self.hbox.Add(self.drawradio, 0, border=3, flag=flags)
        self.hbox.Add(self.drawradio3, 0, border=3, flag=flags)
        self.hbox.Add(self.drawradio2, 0, border=3, flag=flags)
        self.hbox.Add(self.drawbuttonBatch, 0, border=3, flag=flags)
        self.hbox.Add(self.drawbuttonSave, 0, border=3, flag=flags)
        self.vbox.Add(self.hbox, 0, flag=wx.ALIGN_LEFT | wx.TOP)

        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)

    def create_menu(self):
        self.menubar = wx.MenuBar()

        menu_file = wx.Menu()
        m_file = menu_file.Append(-1, "&File \tCtrl-F", "Open a new file")
        self.Bind(wx.EVT_MENU, self.on_file, m_file)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.on_Close, m_exit)

        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&About\tF1", "About the demo")
        # self.Bind(wx.EVT_MENU, self.on_about, m_about)

        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_help, "&Help")
        self.SetMenuBar(self.menubar)

    def on_file(self, cos):
        # ta metoda obsluguje File w Menu
        dlg = wx.FileDialog(self, "Choose a file", self.dirname, "", "*.*", wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            sciezka = os.path.join(dirname, self.filename)
            print(self.filename)
            if self.dirname != dirname:
                self.lista_plikow = ListaDwustronna(self.get_pliki(dirname))
                self.dirname = dirname
                self.zalatw_plik_wynikow(dirname)
                # trzeba tez dopisac zapisywanie przy zmianie katalogu
        dlg.Destroy()
        global ax, ax3
        ax = fig.add_subplot(1, 2, 1)
        ax3 = fig.add_subplot(1, 2, 2)
        self.rysuj_startowa()
        # sprawdzam, ktorym elementem listy plikow jest podana nazwa
        czy_prawidlowy_typ = self.lista_plikow.rewind_to_position(self.filename)
        if czy_prawidlowy_typ:
            self.rysuj(self.filename)
            self.czy_ruszane = True
        else:
            self.rysuj_on_error(self.filename)

    def rysuj(self, nazwa_pliku):
        print(nazwa_pliku)
        nazwa = nazwa_pliku  ## do wykorzystania przy wyswietlaniu na czarno/czerwono
        # nazwa_pliku_sciezka=os.path.join(self.dirname, nazwa_pliku)
        # ------------------------------------------------------------------
        # Decyzja o funkcji znajdujacej diast/syst
        if self.drawradio.GetValue():
            znajdzsyst = self.znajdz_systole_diastole_second_deriv
            self.diastolemethod = "Second Deriv"
        if self.drawradio2.GetValue():
            znajdzsyst = self.znajdz_systole_diastole
            self.diastolemethod = "Sphygmocore"
        if self.drawradio3.GetValue():
            znajdzsyst = self.znajdz_systole_diastole_sphygmocor
            self.diastolemethod = "D Notch"
        # ------------------------------------------------------------------
        central, radial, czas = self.laduj(nazwa_pliku)
        if self.drawbuttonCheckboxNoD.GetValue():
            radial = radial - radial[0]
            print("Bez diastoli")
        ax.clear()
        ax3.clear()
        radial_sys, czas_sys, radial_dia, czas_dia = znajdzsyst(radial, czas)
        TN = czas_dia[0]
        P0 = radial_sys[0]
        # ------------------------------------------------------------------
        # Decyzja czy brac 2/3 czy calosc - uzywam tylko leastsquares i juz
        if self.drawbuttonCheckbox2.GetValue():
            znajdz_PinftauPTN = self.znajdz_PinfRCPTN2_3
            self.falloff = "part"
        else:
            znajdz_PinftauPTN = self.znajdz_PinfRCPTN
            self.falloff = "full"
        PTN, Pinf, b = znajdz_PinftauPTN(radial_dia, czas_dia, TN)
        a = self.znajdz_a(radial_sys, czas_sys, PTN, TN, Pinf, b)
        reservoir_p = self.P_reservoir_fun(radial, czas, a, b, P0, Pinf, TN)
        excess_p = radial - reservoir_p
        excess_pressure = radial - reservoir_p + radial[0]
        czas = 1000 * czas
        wykres = ax.plot(czas, radial, "-k", label="radial pressure")
        ax.plot(czas, reservoir_p, "--k", label="reservoir pressure")
        ax.plot(czas, excess_pressure, "-.k", label="excess pressure")
        if self.show_diastole:
            ax.plot([czas[TN * 128]], radial[TN * 128], ".r", markersize=20)
        ax.set_xlim([0, czas[-1]])
        ax.set_xlabel("time (ms)")
        ax.set_ylabel("pressure (mmHg)")
        ax.legend()

        try:
            czas_fill = arange(0, where(reservoir_p > radial)[0][0])
        except IndexError:
            czas_fill = arange(0, len(czas))
        ax3.plot(czas, central, "-k", label="radial pressure", lw=2)
        ax3.set_xticks([0, 200, 400, 600, 800])
        ax3.set_xlabel("time (ms)")
        ax3.set_ylabel("pressure (mmHg)")

        ## tu zdefiniujemy odstepy, ktore potem beda uzyte na rysunku
        przegiecie = self.znajdz_przegiecie(central)
        y_p = array([min(central), central[przegiecie], max(central)])
        mnoznik = 1000
        odsun = 0.05 * mnoznik
        odstep_strzalka = 0.5
        ujemne = 0.35 * mnoznik
        t_p_min = array([0, 2 * odsun, 0]) - ujemne
        t_p_max = [
            0 - odsun,
            czas[przegiecie] - odsun,
            czas[where(central == max(central))[0][0]] - odsun,
        ]

        ## tutaj segmenty poziome
        ax3.hlines(y_p, t_p_min, t_p_max)
        ax3.set_xlim([-1 * ujemne - odsun, max(czas)])

        ## tutaj strzalki
        ax3.add_patch(
            FancyArrowPatch(
                (-1 * ujemne + odsun, min(central) + odstep_strzalka),
                (-1 * ujemne + odsun, (max(central)) - odstep_strzalka),
                arrowstyle="<->",
                mutation_scale=80,
            )
        )
        ax3.add_patch(
            FancyArrowPatch(
                ((-1 * ujemne + odsun) / 2, min(central) + odstep_strzalka),
                ((-1 * ujemne + odsun) / 2, (central[przegiecie]) - odstep_strzalka),
                arrowstyle="<->",
                mutation_scale=80,
            )
        )
        ax3.add_patch(
            FancyArrowPatch(
                ((-1 * ujemne + odsun) / 2, (central[przegiecie]) + odstep_strzalka),
                ((-1 * ujemne + odsun) / 2, max(central) - odstep_strzalka),
                arrowstyle="<->",
                mutation_scale=80,
            )
        )

        ## tutaj opis strzalek
        ax3.text(
            -1 * ujemne + 1.5 * odsun,
            (min(central) + max(central)) / 2,
            "PP",
            fontsize=15,
        )
        ax3.text(
            (-1 * ujemne + odsun) / 2 + 0.5 * odsun,
            (central[przegiecie] + min(central)) / 2,
            "P1",
            fontsize=15,
        )
        ax3.text(
            (-1 * ujemne + odsun) / 2 + 0.5 * odsun,
            (central[przegiecie] + max(central)) / 2 - 0.5,
            "AP",
            fontsize=15,
        )
        ax3.text(
            czas[where(central == max(central))[0][0]] + odsun,
            max(central) - 1,
            "SBP",
            fontsize=15,
        )
        ax3.text(czas[-1] - 3 * odsun, min(central), "DBP", fontsize=15)
        ### tutaj wypiszemy nazwe pliku - kolor czerwony znaczy, ze plik przetwarzany jest po raz pierwszy, kolor czarny oznacza, ze juz byl obrabiany
        print("rysuje")
        self.canvas.draw()
        ## w tym miejscu wpiszemy wszystkie wyniki do slownika wynikow
        self.results_dictionary[self.filename]["Visited"] = "Yes"
        ## pelna linia radial
        wyniki = self.cechy_linii(radial, czas, TN)
        self.results_dictionary[self.filename]["rad.max"] = str(wyniki[0])
        self.results_dictionary[self.filename]["rad.int"] = str(wyniki[1])
        self.results_dictionary[self.filename]["rad.sys.int"] = str(wyniki[2])
        self.results_dictionary[self.filename]["rad.dia.int"] = str(wyniki[3])
        wyniki = self.cechy_linii(reservoir_p, czas, TN)
        self.results_dictionary[self.filename]["res.max"] = str(wyniki[0])
        self.results_dictionary[self.filename]["res.int"] = str(wyniki[1])
        self.results_dictionary[self.filename]["res.sys.int"] = str(wyniki[2])
        self.results_dictionary[self.filename]["res.dia.int"] = str(wyniki[3])
        wyniki = self.cechy_linii(excess_p, czas, TN)
        self.results_dictionary[self.filename]["ex.max"] = str(wyniki[0])
        self.results_dictionary[self.filename]["ex.int"] = str(wyniki[1])
        self.results_dictionary[self.filename]["ex.sys.int"] = str(wyniki[2])
        self.results_dictionary[self.filename]["ex.dia.int"] = str(wyniki[3])
        #         wyniki=self.cechy_linii(radial-radial[0], czas TN)
        #         self.results_dictionary[self.filename]['rad.pti.eng']=str(wyniki[1])
        #         self.results_dictionary[self.filename]['rad.sys.pti.eng']=str(wyniki[2])
        #         self.results_dictionary[self.filename]['rad.dia.pti.eng']=str(wyniki[3])
        #         wyniki=self.cechy_linii(reservoir_p-radial[0], czas TN)
        #         self.results_dictionary[self.filename]['res.pti.eng']=str(wyniki[1])
        #         self.results_dictionary[self.filename]['res.sys.pti.eng']=str(wyniki[2])
        #         self.results_dictionary[self.filename]['res.dia.pti.eng']=str(wyniki[3])
        # TUTU
        self.results_dictionary[self.filename]["diast"] = self.diastolemethod
        self.results_dictionary[self.filename]["falloff"] = self.falloff

    def rysuj_batch(self, nazwa_pliku):
        print(nazwa_pliku)
        nazwa = nazwa_pliku  ## do wykorzystania przy wyswietlaniu na czarno/czerwono
        # nazwa_pliku_sciezka=os.path.join(self.dirname, nazwa_pliku)
        # ------------------------------------------------------------------
        # Decyzja o funkcji znajdujacej diast/syst
        if self.drawradio.GetValue():
            znajdzsyst = self.znajdz_systole_diastole_second_deriv
        if self.drawradio2.GetValue():
            znajdzsyst = self.znajdz_systole_diastole
        if self.drawradio3.GetValue():
            znajdzsyst = self.znajdz_systole_diastole_sphygmocor
        # ------------------------------------------------------------------
        central, radial, czas = self.laduj(nazwa_pliku)
        if self.drawbuttonCheckboxNoD.GetValue():
            radial = radial - radial[0]
            print("Bez diastoli")
        # ax.clear()
        # ax2.clear()
        radial_sys, czas_sys, radial_dia, czas_dia = znajdzsyst(radial, czas)
        TN = czas_dia[0]
        P0 = radial_sys[0]
        # ------------------------------------------------------------------
        # Decyzja czy brac 2/3 czy calosc - uzywam tylko leastsquares i juz
        if self.drawbuttonCheckbox2.GetValue():
            znajdz_PinftauPTN = self.znajdz_PinfRCPTN2_3
        else:
            znajdz_PinftauPTN = self.znajdz_PinfRCPTN
        PTN, Pinf, b = znajdz_PinftauPTN(radial_dia, czas_dia, TN)
        a = self.znajdz_a(radial_sys, czas_sys, PTN, TN, Pinf, b)
        reservoir_p = self.P_reservoir_fun(radial, czas, a, b, P0, Pinf, TN)
        excess_p = radial - reservoir_p
        excess_pressure = radial - reservoir_p + radial[0]
        ## w tym miejscu wpiszemy wszystkie wyniki do slownika wynikow
        self.results_dictionary[self.filename]["Visited"] = "Yes"
        ## pelna linia radial
        wyniki = self.cechy_linii(radial, czas, TN)
        self.results_dictionary[self.filename]["rad.max"] = str(wyniki[0])
        self.results_dictionary[self.filename]["rad.int"] = str(wyniki[1])
        self.results_dictionary[self.filename]["rad.sys.int"] = str(wyniki[2])
        self.results_dictionary[self.filename]["rad.dia.int"] = str(wyniki[3])
        wyniki = self.cechy_linii(reservoir_p, czas, TN)
        self.results_dictionary[self.filename]["res.max"] = str(wyniki[0])
        self.results_dictionary[self.filename]["res.int"] = str(wyniki[1])
        self.results_dictionary[self.filename]["res.sys.int"] = str(wyniki[2])
        self.results_dictionary[self.filename]["res.dia.int"] = str(wyniki[3])
        wyniki = self.cechy_linii(excess_p, czas, TN)
        self.results_dictionary[self.filename]["ex.max"] = str(wyniki[0])
        self.results_dictionary[self.filename]["ex.int"] = str(wyniki[1])
        self.results_dictionary[self.filename]["ex.sys.int"] = str(wyniki[2])
        self.results_dictionary[self.filename]["ex.dia.int"] = str(wyniki[3])
        self.results_dictionary[self.filename]["diast"] = self.diastolemethod
        self.results_dictionary[self.filename]["falloff"] = self.falloff

    def rysuj_startowa(self):
        ax.plot([0, 1], [0, 1], color="white")
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        ax.text(0.3, 0.7, "This is the WCT-AVA program", fontsize=22)
        ax.text(
            0.28, 0.6, "Go to the File menu and select a file to start", fontsize=16
        )
        ax.text(0.44, 0.5, "version 0.01", fontsize=16)
        ax.text(0.3, 0.4, "Jaroslaw Piskorski, jaropis@zg.home.pl", fontsize=16)

    def rysuj_on_error(self):
        ax.clear()
        ax.plot([0, 1], [0, 1], color="white")
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        ax.text(0.35, 0.7, "WRONG FILE TYPE!!!", fontsize=22)
        ax.text(0.28, 0.6, "Go to the File menu and select another file", fontsize=16)
        self.canvas.draw()

    def rysuj_batch_info(self):
        ax.clear()
        ax.plot([0, 1], [0, 1], color="white")
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        ax.text(0.0, 0.7, "batch mode", fontsize=22)
        # self.canvas.draw()
        ax2.clear()
        ax2.plot([0, 1], [0, 1], color="white")
        ax2.get_xaxis().set_ticks([])
        ax2.get_yaxis().set_ticks([])
        ax2.text(0.0, 0.7, "batch mode", fontsize=22)
        self.canvas.draw()

    def on_Next(self, cos):
        self.filename = self.lista_plikow.get_next()
        self.rysuj(self.filename)

    def on_Previous(self, cos):
        self.filename = self.lista_plikow.get_previous()
        self.rysuj(self.filename)

    def on_Repeat(self, cos):
        self.rysuj(self.filename)

    def on_CheckBox(self, cos):
        sender = cos.GetEventObject()
        isChecked = sender.GetValue()
        if isChecked:
            self.show_diastole = True
        else:
            self.show_diastole = False
        print(self.show_diastole, "checkboxxxxxx")

    def on_Batch(self, cos):
        lista_plikow = self.get_pliki(self.dirname)
        self.rysuj_batch_info()
        for plik in lista_plikow:
            self.filename = plik
            self.rysuj_batch(plik)
        self.rysuj(self.filename)

    def on_Save(self, cos):
        nazwa_obrazka = self.filename[0:-4] + ".png"
        zapisz_tak = os.path.join(self.dirname, nazwa_obrazka)
        savefig(zapisz_tak, dpi=600, format="png")

    def on_Close(self, cos):
        if self.czy_ruszane:
            plik = open("1results.csv", "w")
            self.create_content_file(plik)
            plik.close()
        else:
            print("nieruszane!!!")
        print("patrzcie!!! KONIEC!!!")
        self.Destroy()
        raise SystemExit

    def laduj(self, nazwa_pliku):
        plik = open(nazwa_pliku)
        plik.readline()
        plik.readline()
        plik.readline()
        plik.readline()
        radial = []
        central = []
        for linijka in plik:
            linijka = linijka.replace(",", ".")
            lin_dane = linijka.split("\t")
            lin_dane = [float(i) for i in lin_dane]
            try:
                radial.append(lin_dane[2])
                central.append(lin_dane[3])
            except IndexError:
                break
            czas = range(0, len(radial))
        central = array(central)
        radial = array(radial)
        czas = array(czas) / 128.0
        plik.close()
        return central, radial, czas

    def znajdz_systole_diastole_sphygmocor(self, radial, czas):
        plik = open(self.filename)
        plik.readline()
        ED = int(plik.readline().split()[-9])  ## ED to pozycja poczatku rozkurczu
        ED = int(round(128 * ED / 1000.0 + 1))
        print(ED, len(radial))
        plik.close()
        radial_sys = radial[0:ED]
        radial_dia = radial[ED:]
        czas_sys = czas[0:ED]
        czas_dia = czas[ED:]
        print("sphygmocor")
        return (
            radial_sys,
            czas_sys,
            radial_dia,
            czas_dia,
        )  ## zwracamy cisnienie na systoli, czas systoli, cisnienie na diastoli, czas diastoli

    def znajdz_systole_diastole(self, radial, czas):
        dlugosc_radiala = len(radial)
        rozniczka = diff(radial)
        ostatni = len(rozniczka)
        licznik_od_konca = 1
        while rozniczka[ostatni - 1] > 0:
            licznik_od_konca += 1
            ostatni -= 1
        licznik_do_przegiecia = 0
        # print "pierwsza petla"
        while rozniczka[ostatni - 1] <= 0:
            licznik_do_przegiecia += 1
            ostatni -= 1
        # print "druga petla"
        while rozniczka[ostatni - 1] > 0:
            licznik_do_przegiecia += 1
            ostatni -= 1
        print("dichrotic notch")
        return (
            radial[0 : dlugosc_radiala - licznik_od_konca - licznik_do_przegiecia + 1],
            czas[0 : dlugosc_radiala - licznik_od_konca - licznik_do_przegiecia + 1],
            radial[dlugosc_radiala - licznik_od_konca - licznik_do_przegiecia :],
            czas[dlugosc_radiala - licznik_od_konca - licznik_do_przegiecia :],
        )  ## zwracamy cisnienie na systoli, czas systoli, cisnienie na diastoli, czas diastoli

    def znajdz_systole_diastole_second_deriv(
        self, radial, czas
    ):  ## wzgledem drugiej pochodnej, uzywajac wyniku ze sphygmokoru, zeby znalezc koniec systolii - podem wybieramy ostatnie zero drugiej pochodnej systolii
        (
            radial_sys_local,
            czas_sys_local,
            radial_dia_local,
            czas_dia_local,
        ) = self.znajdz_systole_diastole_sphygmocor(radial, czas)
        A = diff(diff(radial_sys_local))
        indeksy_zerowe = where(A[:-1] * A[1:] < 0)
        TN2 = indeksy_zerowe[0][-1] + 1
        return (
            radial[0:TN2],
            czas[0:TN2],
            radial[TN2:],
            czas[TN2:],
        )  ## zwracamy cisnienie na systoli, czas systoli, cisnienie na diastoli, czas diastoli

    # def znajdz_systole_diastole_second_deriv(self, radial, czas):## wzgledem drugiej pochodnej
    #     dlugosc_radiala=len(radial)
    #     TN2=where(diff(radial)==min(diff(radial)))[0] ## tutaj mamy druga definicje czasu poczatku rozkurczu
    #     print "second deriv"
    #     return radial[0:TN2], czas[0:TN2], radial[TN2:], czas[TN2:]## zwracamy cisnienie na systoli, czas systoli, cisnienie na diastoli, czas diastoli

    def znajdz_PinfRCPTN(self, radial_dia, czas_dia, TN):  ## leastsquare
        funkcja = (
            lambda p: (p[0] - p[1]) * exp(-(czas_dia - TN) * p[2]) + p[1]
        )  # p[0] - PTN, p[1] - Pinf, p[2] - b
        print(len(czas_dia))
        error_fun = lambda p, radial_dia: funkcja(p) - radial_dia
        starting = [radial_dia[0], 0.0, 1.0]
        p, success = leastsq(error_fun, starting, args=(radial_dia), maxfev=10000)
        print(success, "RCPTN leastsquares dla full fall-off")
        return p[0], p[1], p[2]

    def znajdz_PinfRCPTN2(self, radial_dia, czas_dia, TN):  ## nedler
        funkcja = (
            lambda p: (p[0] - p[1]) * exp(-(czas_dia - TN) * p[2]) + p[1]
        )  # p[0] - PTN, p[1] - Pinf, p[2] - b
        error_fun = lambda p: sum((funkcja(p) - radial_dia) ** 2)
        starting = [radial_dia[0], 0.0, 1.0]
        p = fmin(
            error_fun, starting, maxfun=1000, xtol=1e-11
        )  # , args=(czas, radial_dia))
        return p

    def znajdz_PinfRCPTN2_3(self, radial_dia, czas_dia, TN):  ## leastsquare
        ## dwie kolejne linijki obliczaja ostatnie 2/3 rozkurczu
        radial_dia = radial_dia[int(round(1 / 3.0 * len(radial_dia))) :]
        czas_dia = czas_dia[int(round(1 / 3.0 * len(czas_dia))) :]
        print(len(czas_dia))
        funkcja = (
            lambda p: (p[0] - p[1]) * exp(-(czas_dia - TN) * p[2]) + p[1]
        )  # p[0] - PTN, p[1] - Pinf, p[2] - b
        error_fun = lambda p, radial_dia: funkcja(p) - radial_dia
        starting = [radial_dia[0], 0.0, 1.0]
        p, success = leastsq(error_fun, starting, args=(radial_dia))
        print(success, "RCPTN leastsquares dla 2/3")
        return p[0], p[1], p[2]

    def znajdz_PinfRCPTN22_3(self, radial_dia, czas_dia, TN):  ## nedler
        funkcja = (
            lambda p: (p[0] - p[1]) * exp(-(czas_dia - TN) * p[2]) + p[1]
        )  # p[0] - PTN, p[1] - Pinf, p[2] - b
        error_fun = lambda p: sum((funkcja(p) - radial_dia) ** 2)
        starting = [radial_dia[0], 0.0, 1.0]
        p = fmin(
            error_fun, starting, maxfun=1000, xtol=1e-11
        )  # , args=(czas, radial_dia))
        return p

    def calka(self, P, a, b, t):
        tau = arange(1 / 128.0, t + 1 / 128.0, 1 / 128.0)
        # print tau
        integral = trapz(a * P * exp((a + b) * tau), dx=1 / 128.0)
        return integral

    def calka2(self, P, a, b, czas):
        integral = simps(a * P * exp((a + b) * czas), czas)
        return integral

    def znajdz_a(self, radial_sys, czas_sys, PTN, TN, Pinf, b):
        P0 = radial_sys[0]  ## cisnienie na poczatku systonle (P(0))
        error_fun = (
            lambda a: PTN
            - b / (a + b) * Pinf
            - exp(-(a + b) * TN)
            * (self.calka2(radial_sys, a, b, czas_sys) + P0 - b * Pinf / (a + b))
        )
        starting = 0.1
        p = fmin(error_fun, starting)
        a = p
        print(
            PTN,
            b / (a + b) * Pinf
            + exp(-(a + b) * TN)
            * (self.calka2(radial_sys, a, b, czas_sys) + P0 - b * Pinf / (a + b)),
        )
        return p

    def P_reservoir_fun(self, radial_sys, czas_sys, a, b, P0, Pinf, TN):
        P_res_sys = [radial_sys[0]]
        for t in czas_sys[1:]:
            P_res_sys.append(
                (
                    b / (a + b) * Pinf
                    + exp(-(a + b) * t)
                    * (
                        self.calka2(
                            radial_sys[0 : int(t * 128) + 1],
                            a,
                            b,
                            czas_sys[0 : int(t * 128) + 1],
                        )
                        + P0
                        - b * Pinf / (a + b)
                    )
                )[0]
            )
        # print P_res_sys[-1]
        return array(P_res_sys)

    def P_reservoir_diastole(radial_dia, czas_dia, b, Pinf, PTN, TN):
        P_res_dia = []
        for t in czas_dia:
            P_res_dia.append(Pinf + (PTN - Pinf) * exp(-b * (t - TN)))
        return array(P_res_dia)

    ## tu bedziemy zbierac cechy linii
    def cechy_linii(self, cisnienie, czas, TN):
        TNinteger = int(TN * 128)
        maks_linii = max(cisnienie)
        calka_pelna = simps(cisnienie, czas)
        calka_sys = simps(cisnienie[0:TNinteger], czas[0:TNinteger])
        calka_dia = simps(cisnienie[TNinteger:], czas[TNinteger:])
        return maks_linii, calka_pelna, calka_sys, calka_dia

    def znajdz_przegiecie(self, central):
        A = diff(diff(central))
        indeksy_zerowe = where(A[:-1] * A[1:] < 0)
        return indeksy_zerowe[0][1]

    def zalatw_plik_wynikow(self, dirname):
        # Ta metoda otwiera plik wynikowy i wczytuje jego zawartosc
        # do slownika self.wyniki. Jezeli plik nie istnieje, to jest tworzony
        # wraz ze struktura, przy czym jest oczywiscie pusty
        try:
            self.result_file.close()  # czy jakis plik wczesniej byl otwierany?
        except AttributeError:
            try:  # czy mozna wczytac istniejacy w katalogu plik?
                self.result_file = open(os.path.join(dirname, "1results.csv"), "r")
            except IOError:  # jezeli plik nie istnieje, to go stworz
                self.result_file = open(os.path.join(dirname, "1results.csv"), "w")
                self.create_results_dictionary()
                self.create_content_file(self.result_file)
                self.result_file.close()  # stworzenie pliku, jezeli nie istnial
                self.result_file = open(os.path.join(dirname, "1results.csv"), "r")
        else:  # co zrobic, jezeli jakis plik byl wczesniej otwierany
            try:  # czy mozna wczytac istniejacy w nowym katalogu plik?
                self.result_file = open(os.path.join(dirname, "1results.csv"), "r")
            except IOError:  # jezeli plik nie istnieje, to go stworz
                self.result_file = open(os.path.join(dirname, "1results.csv"), "w")
                self.create_results_dictionary(self.lista_plikow.lista)
                self.create_content_file(self.result_file)
                self.result_file.close()  # stworzenie pliku, jezeli nie istnial
                self.result_file = open(os.path.join(dirname, "1results.csv"), "r")
        finally:
            print(self.result_file)
            self.get_content_file_results(
                self.result_file
            )  # TUTU musisz oprogramowac te linijke!!!
            print("odczytalem")
            self.result_file.close()

    def get_content_file_results(self, plik):
        dane = reader(plik, delimiter=";")
        dane = list(dane)
        print(type(self.results_dictionary))
        for i in range(1, len(dane)):
            self.results_dictionary[dane[i][0]] = {
                "Visited": dane[i][1],
                "rad.max": dane[i][2],
                "rad.int": dane[i][3],
                "rad.sys.int": dane[i][4],
                "rad.dia.int": dane[i][5],
                "res.max": dane[i][6],
                "res.int": dane[i][7],
                "res.sys.int": dane[i][8],
                "res.dia.int": dane[i][9],
                "ex.max": dane[i][10],
                "ex.int": dane[i][11],
                "ex.sys.int": dane[i][12],
                "ex.dia.int": dane[i][13],
                "diast": self.diastolemethod,
                "falloff": self.falloff,
            }

    def create_results_dictionary(self):
        for nagranie in self.lista_plikow.lista:
            self.results_dictionary[nagranie] = {
                "Visited": "No",
                "rad.max": "NA",
                "rad.int": "NA",
                "rad.sys.int": "NA",
                "rad.dia.int": "NA",
                "res.max": "NA",
                "res.int": "NA",
                "res.sys.int": "NA",
                "res.dia.int": "NA",
                "ex.max": "NA",
                "ex.int": "NA",
                "ex.sys.int": "NA",
                "ex.dia.int": "NA",
                "diast": "NA",
                "falloff": "NA",
            }

    def create_content_file(self, plik):
        plik.write(
            "Filename;Visited;rad.max;rad.int;rad.sys.int;rad.dia.int;res.max;res.int;res.sys.int;res.dia.int;ex.max;ex.int;ex.sys.int;ex.dia.int;diast;falloff\n"
        )
        for nazwa in sorted(self.results_dictionary.keys()):
            lancuch_wynikowy = (
                nazwa
                + ";%(Visited)s;%(rad.max)s;%(rad.int)s;%(rad.sys.int)s;%(rad.dia.int)s;%(res.max)s;%(res.int)s;%(res.sys.int)s;%(res.dia.int)s;%(ex.max)s;%(ex.int)s;%(ex.sys.int)s;%(ex.dia.int)s;%(diast)s;%(falloff)s\n"
                % self.results_dictionary[nazwa]
            )
            plik.write(lancuch_wynikowy)


if __name__ == "__main__":
    app = wx.PySimpleApp()
    app.frame = ReservoirFrame()
    app.frame.Show()
    app.MainLoop()
