'''
For the gene info, it might be a good time to introduce SQL, the relational database manager. Basically we need one table A 
that stores gene -> PAM, and another table B that stores PAM -> variants and their scores. This will make the counting easier 
by only looking at table A. I haven’t thought enough about what would be the best database layout for croton, but we are going 
to use SQL anyways so might be better to start looking at how SQL works with Python
    - https://docs.python.org/3/library/sqlite3.html
'''

# A gene PAM [coord of PAM]
# PAM variant [coordinate of variant, id, column whether disrupt]
# variant prediction [for all stats]
# Reference 

#N coord => PAM
#ref allele, alt allele, allele 


import sqlite3

def make_PAMdb():
    conn = sqlite3.connect('pam.db')
    c = conn.cursor()
    # Create table
    c.execute('''CREATE TABLE pam
                (genename, pam)''')

    # Insert a row of data
    c.execute("INSERT INTO pam VALUES ('ABCCD', 'fwejfijw')")

    # Save (commit) the changes
    conn.commit()

    # We can also close the connection if we are done with it.
    # Just be sure any changes have been committed or they will be lost.
    conn.close()

def make_vardb():
    conn = sqlite3.connect('var.db')
    c = conn.cursor()
    # Create table
    c.execute('''CREATE TABLE variant
                (genename, var)''')

    # Insert a row of data
    c.execute("INSERT INTO variant VALUES ('ABCCD', 'fwejfijw')")

    # Save (commit) the changes
    conn.commit()

    # We can also close the connection if we are done with it.
    # Just be sure any changes have been committed or they will be lost.
    conn.close()

#Make tuples
'''
# Never do this -- insecure!
symbol = 'RHAT'
c.execute("SELECT * FROM stocks WHERE symbol = '%s'" % symbol)

# Do this instead
t = ('RHAT',)
c.execute('SELECT * FROM stocks WHERE symbol=?', t)
print(c.fetchone())

# Larger example that inserts many records at a time
purchases = [('2006-03-28', 'BUY', 'IBM', 1000, 45.00),
             ('2006-04-05', 'BUY', 'MSFT', 1000, 72.00),
             ('2006-04-06', 'SELL', 'IBM', 500, 53.00),
            ]
c.executemany('INSERT INTO stocks VALUES (?,?,?,?,?)', purchases)
'''