# Search Disease App

## Requirements

- Python 3.8 (3.8.16)

## Local development

- Add secrets file `.streamlit/secrets.toml`
```toml
[app_credentials]
name = "Example Example"
password = "example"
username = "root"
```

- To Run project

```shell
streamlit run .\app.py
```

## Reference

- Style for my web app: https://share.streamlit.io/streamlit/example-app-bert-keyword-extractor/main/app.py

## Utils

- https://fsymbols.com/generators/tarty/

## Docs

- ขั้นตอน data preparation
- (js)รัน -> output -> เอาไปใช้ยังไงต่อ (block diagram)
- jupyternotebook -> detail แต่ละ data source

- เขียนขั้นตอนการ deploy
- flow ของการดึงข้อมูล - ของแต่ละตารางดึงข้อมูลมาอย่างไร
- อธิบายแต่ละ part ของ output ที่เอามาใช้งาน (วิธีใช้ไม่ต้อง)
- \*จาก issue เรื่อง data dup -> ให้ uppercase ฝั่งขาที่จะเอามา search
