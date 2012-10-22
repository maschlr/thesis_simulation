# Script to read the weater data supported the Universidad de Piura in some ugly format
import sqlite3
import csv
import os
import re
from numpy import array, floor

class weatherData(object):
    def __init__(self):
        self.data = []
        self.listOfFiles=[]
        self.year = None

    def setYear(self, year):
        self.year=year

    def getAllWeatherFiles(self,directory=os.path.join(os.getcwd(),'weatherfiles/')):
        dataFiles = os.listdir(directory)
        dataFiles.sort()
        for file in dataFiles:
            print file
        answer = raw_input('Should these files be scanned for weather data? (y/n)')
        if answer == 'no' or answer == 'n':
            answer = raw_input('Do you want to pop some out of the list or scan another directory? (p/a)')
            if answer == 'p' or answer == 'pop':
                for i in range(len(dataFiles)):
                    print i, ': ', dataFiles[i]
                pop = raw_input('Please type in the index of the files that should be excluded, separated by commata: ')
                n = 0
                for i in pop.split(','):
                    dataFiles.pop(int(i)-n)
                    n+=1
            elif answer == 'a':
                directory = raw_input('Which directory should be scanned for files instead? \n Full path: ')
                self.getAllWeatherFiles(directory)
        elif answer=='yes' or answer=='y':
            pass
        else: 
            print("You didn't answer the question properly, we'll try another time!")
            self.getAllWeatherFiles(directory)
        self.listOfFiles = []
        for dataFile in dataFiles:
            self.listOfFiles.append(os.path.join(directory, dataFile))

    def readDataFromFiles(self):
        self.data = []
        for file in self.listOfFiles:
            f = open(file, 'r')
            f.readline()
            f.readline()
            yearMinutes = 0
            i = 1
            for l in f.readlines():
                l = l.split()
                l[0] = l[0].split('/')
                l[1] = l[1].split(':')
                if i < 2:
                    yearMinutes = self.getYearMinutes([int(b) for b in l[0]])
                    yearMinutes += int(l[1][0])*60
                    self.setYear(int(l[0][2]))
                # The colums are:(day, month, year, hours, minutes, yearMinutes, temperature, humidity, windspeed, insolation) 
                self.data.append([int(l[0][0]),int(l[0][1]),int(l[0][2]),int(l[1][0]),int(l[1][1]),yearMinutes+i*30,float(l[2]),float(l[5])/100,float(l[7]),float(l[19])])
                i+=1
            f.close()
    
    def getYearMinutes(self, dayMonthYear):
        #Computes the minutes passed since 01.Jan 0:00
        hoursTotal = 0
        leapyears = range(1900,2200,4)
        leapyears.remove(1900)
        leapyears.remove(2100)
        daysInMonth = [31, 28, 31, 30, 31, 30, 31 ,31, 30, 31, 30, 31]
        if dayMonthYear[2] in leapyears:
            daysInMonth[1] = 29
        if dayMonthYear[1]>1:
            for i in range(dayMonthYear[1]-1):
                hoursTotal += daysInMonth[i]*24
        return (hoursTotal+(dayMonthYear[0]-1)*24)*60

    def writeFile(self, filename='default.csv', path=''):
        f = open(os.path.join(path,filename), 'w')
        header = ('day', 'month', 'year', 'hours', 'minutes', 'yearMinutes', 'temperature', 'humidity', 'windspeed', 'insolation')
        csv_writer = csv.writer(self.f)
        csv_writer.writerow(header)
        for row in self.data:
            csv_writer.writerow(row)
        f.close()

    def writeDB(self, dbfile='weather.db', path='', tableName='weatherData'):
        conn = sqlite3.connect(os.path.join(path, dbfile))
        db = conn.cursor()

        db.execute("CREATE TABLE IF NOT EXISTS %s (day int, month int, year int, hours int, \
                minutes int, yearMinutes int, temperature real, humidity real, windspeed real, insolation real)" %tableName)
        for row in self.data:
            # create a tuple to pass it as argument to the execute function
            db.execute('INSERT INTO %s VALUES (?,?,?,?,?,?,?,?,?,?)' %tableName, tuple(row)) 
        conn.commit()
        db.close()

    def readDB(self, month, startDay, startTime, timeSpan, year=None, dbfile='weather.db', path='', tableName='weatherData'):
        startMinutes = self.getYearMinutes([startDay, month, (year or self.year)])+startTime*60
        endMinutes = startMinutes+timeSpan*60
        conn = sqlite3.connect(os.path.join(path,dbfile))
        db = conn.cursor()
        db.execute('SELECT * FROM %s WHERE month=? AND yearMinutes>=? AND yearMinutes<=? ORDER BY yearMinutes' %tableName, (month, startMinutes, endMinutes) )
        self.relevantData = db.fetchall()
        db.close()

        


if __name__ == '__main__':
    print('Welcome to the data reading script that takes the ugly files,') 
    print('provided by the Piura weather radar station and gives you    ')
    print('either a .csv file or a sqlite3 database with the numbers that')
    print('really matter without all this useless shit. ')
    print('Your output format do you want to use?')
    inputChoice = raw_input(' 1 : .csv   2 : sqlite3   3 : both ')
    print('')
    print('Do you want to search the current directory for files?')
    directorySearch = raw_input('[y/n] ')
    # Only get the files with the known pattern
    dataFiles = filter(lambda x: re.match('[0-9]{6}\.txt', x), os.listdir(os.getcwd()))
    dataFiles.sort()
    
    print('')
    print('I have found the following files with potential data:')
    print(dataFiles)
    switchContinue = raw_input('Continue? [y/n] ')
    wd = weatherData(dataFiles)
    wd.writeFile()
    wd.writeDB()




