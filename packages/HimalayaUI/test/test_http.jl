using Test, HTTP, SQLite, JSON3
using HimalayaUI

function with_test_server(f, db::SQLite.DB)
    port = HimalayaUI.find_free_port()
    HimalayaUI.start_test_server!(db, port)
    try
        f(port, "http://127.0.0.1:$port")
    finally
        HimalayaUI.stop_test_server!()
    end
end
