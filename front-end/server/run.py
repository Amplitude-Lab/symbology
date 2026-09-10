from __future__ import annotations

import os
import uvicorn

from app.config import HOST, PORT

if __name__ == "__main__":
    uvicorn.run("app.main:app", host=HOST, port=PORT, log_level="info", reload=os.environ.get("SYMBOLOGY_RELOAD") == "1", reload_dirs=["app"])
