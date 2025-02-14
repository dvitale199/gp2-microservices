from fastapi import FastAPI, HTTPException
from carriers import process_carriers_data
from models import CarriersResponse

app = FastAPI(
    title="GP2 Carriers API",
    description="API for processing GP2 carriers data",
    version="1.0.0"
)

@app.post("/process-carriers", response_model=CarriersResponse)
async def process_carriers():
    try:
        result = process_carriers_data()
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))