#ifndef SKYLENS_SQLITEDB_H
#define SKYLENS_SQLITEDB_H

#include <string>
#include <sqlite3.h>

namespace skylens {  

  /// DB accessor class for SQLite.
  class SQLiteDB {
  public:
    /// Constructor.
    SQLiteDB();
    /// Destructor.
    ~SQLiteDB();
    /// Result class for SQLiteDB.
    class SQLiteResult {
    public:
      /// Constructor.
      SQLiteResult();
      /// Destructor.
      ~SQLiteResult();
      /// Get number of rows in result set.
      unsigned int getRowCount();
      /// Get number of fields in each row of the result set.
      unsigned int getFieldCount();
      /// Get the name of the field \p i.
      char* getFieldName(unsigned int i);
      /// Access the complete row of the result set.
      /// For result sets with multiple rows, call this function
      /// until it returns \p NULL.
      char** getRow();
      friend class SQLiteDB;
    private:
      char** table;
      int nrows, ncols, counter;
    };

    /// Connect to DB \p dbfile.
    /// If \p dbfile is <tt>:memory:</tt>, then a private, 
    /// temporary in-memory database is created for the connection.
    void connect(std::string dbfile);
    /// Execute a query using \p callback.
    /// See http://www.sqlite.org/c3ref/exec.html and
    /// http://www.sqlite.org/quickstart.html for details.
    void exec(std::string query, int (*callback)(void*,int,char**,char**));
    /// Retrieve data from query.
    SQLiteResult query(std::string query);
    /// Type of DBHandle (returns 0).
    int getDBType();
    /// Direct access to the DB connection.
    sqlite3* db;
    void checkRC(int rc);
  private:
    char *zErrMsg;
    std::string database;
    
  };

} // end namespace

#endif
