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

    def getAllWeatherFiles(self,directory=os.path.join(os.getcwd(),'weatherfiles/')):
        dataFiles = filter(lambda x: re.match('[0-9]{6}\.txt', x), os.listdir(directory))
        dataFiles.sort()
        self.listOfFiles = []
        for dataFile in dataFiles:
            self.listOfFiles.append(os.path.join(directory, dataFile))

    def readDataFromFiles(self):
        self.data = []
        for file in self.listOfFiles:
            f = open(file, 'r')
            f.readline()
            f.readline()
            i = 1
            for l in f.readlines():
                l = l.split()
                l[0] = l[0].split('/')
                l[1] = l[1].split(':')
                # The colums are:(day, month, year, hours, minutes, seconds, temperature, humidity, windspeed, insolation) 
                self.data.append([int(l[0][0]),int(l[0][1]),int(l[0][2]),int(l[1][0]),int(l[1][1]),i*1800,float(l[2]),float(l[5])/100,float(l[7]),float(l[19])])
                i+=1
            f.close()

    def writeFile(self, filename='default.csv', path=''):
        f = open(os.path.join(path,filename), 'w')
        header = ('day', 'month', 'year', 'hours', 'minutes', 'seconds', 'temperature', 'humidity', 'windspeed', 'insolation')
        csv_writer = csv.writer(self.f)
        csv_writer.writerow(header)
        for row in self.data:
            csv_writer.writerow(row)
        f.close()

    def writeDB(self, dbfile='default.db', path='', tableName='weatherData'):
        conn = sqlite3.connect(os.path.join(path, dbfile))
        db = conn.cursor()

        db.execute("CREATE TABLE IF NOT EXISTS %s (day int, month int, year int, hours int, \
                minutes int, seconds int, temperature real, humidity real, windspeed real, insolation real)" %tableName)
        for row in self.data:
            # create a tuple to pass it as argument to the execute function
            db.execute('INSERT INTO %s VALUES (?,?,?,?,?,?,?,?,?,?)' %tableName, tuple(row)) 
        conn.commit()
        db.close()

    def readDB(self, month, startDay, startTime, timeSpan, dbfile='default.db', path='', tableName='weatherData'):
        startSeconds=(startDay*24+startTime)*3600
        endSeconds = startSeconds+timeSpan*3600
        conn = sqlite3.connect(os.path.join(path,dbfile))
        db = conn.cursor()
        db.execute('SELECT * FROM %s WHERE month=? AND seconds>=? AND seconds<=? ORDER BY seconds' %tableName, (month, startSeconds, endSeconds) )
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




