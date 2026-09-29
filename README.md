# FME CSMap Pipeline

DEM（数値標高モデル）からCS立体図（長野県林業総合センター考案の可視化手法）を生成する
スタンドアロンCLIツールです。FMEのSystemCaller等から呼び出す用途を想定しています。

本体は [`fme_csmap_pipeline/`](fme_csmap_pipeline/) フォルダーです。

詳細な機能説明・入力形式・設定項目・既知の制約は以下を参照してください。

- [fme_csmap_pipeline/README_ja.md](fme_csmap_pipeline/README_ja.md) — 機能・使い方の詳細
- [fme_csmap_pipeline/INPUT_GUIDE_ja.md](fme_csmap_pipeline/INPUT_GUIDE_ja.md) — 対応入力形式ごとの詳細
- [fme_csmap_pipeline/CHANGELOG.md](fme_csmap_pipeline/CHANGELOG.md) — バージョンごとの変更点
- [fme_csmap_pipeline/VALIDATION.txt](fme_csmap_pipeline/VALIDATION.txt) — 検証記録・未確認事項

姉妹プロジェクトとして、同じ処理をQGIS上で行うプラグイン
[csmap-sheets](https://github.com/yo-1/csmap-sheets) があります。
両者はコード・バージョンとも独立して保守されており、機能が完全に同期しているとは限りません。

## ライセンス

GNU GPL v3.0 only。詳細は [LICENSE](LICENSE) を参照してください。

Copyright (C) 2026 Yoichi Wada.
