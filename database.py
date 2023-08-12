import sqlite3

# Connect to the SQLite database
conn = sqlite3.connect('Zombies\db.sqlite3')
c = conn.cursor()

c.execute('CREATE TABLE IF NOT EXISTS synonyms (drugID INTEGER, synonym TEXT, '
               'FOREIGN KEY (drugID) REFERENCES drugsList(id))')
conn.commit()